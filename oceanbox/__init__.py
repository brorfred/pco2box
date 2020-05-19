import time
import copy

import numpy as np
import pylab as pl
import pandas as pd

import seawater as sw
from .wind2pv import wind2pv, schmidt
from . import stload
from . import reuerflux
import co2sys

try:
    import click
    HAS_CLICK = True
except ImportError:
    HAS_CLICK = False
    
def nodim(a,b):  return (a-b) / ((a+b)/2)
def nodimb(a,b): return (a-b) / b

def clean(a):
    a[np.isinf(a)]=np.nan
    return a[~isnan(a)]


class BoxModel(object):
    """1D boxmodel to calculate air-sea fluxes of O2 and CO2"""
    def __init__(self, svec=365):
        self.zerovec = np.squeeze(np.zeros(svec))
        self.tvec = np.arange(len(self.zerovec[:,...]))
        self.reset()

    def reset(self):
        self.setup_physics()
        self.setup_o2ar()
        self.setup_carb()
        
    def setup_physics(self):
        initpar = {"temp":28.0, "salt":33.0, "nwnd":10.0, "mldp":30.0,
                   "ncpm":0.0,  "ncpc":0.0,  "bubl":0.01}
        for key in initpar:
            setattr(self, key, self.zerovec + initpar[key])

    def setup_o2ar(self):
        """Setup needed parameters to model O2 and Ar"""
        for key in ['o2armsk', 'o2fl', 'arfl']:
            setattr(self, key, self.zerovec.copy())
        self.o2ct = sw.satO2(self.salt, self.temp)
        self.arct = sw.satAr(self.salt, self.temp)
        self.o2st = sw.satO2(self.salt, self.temp)
        self.arst = sw.satAr(self.salt, self.temp)
        #self.o2pv =  wind2pv(self.temp, self.nwnd)
    
    def setup_carb(self, pco2_atm=None, pco2=None, talk=None):
        """Setup needed parameters to model CO2"""
        self.pco2_atm = self.zerovec+395.0 if pco2_atm is None else pco2_atm
        self.pco2     = self.zerovec+395.0 if pco2     is None else pco2.copy()
        if talk is None:
            self.talk = self.zerovec + co2sys.calc_TAlk(self.salt, self.temp)
        else:
            self.talk = talk
        co = co2sys.CarbonateSystem(
            self.salt, self.temp, TA=self.talk, pCO2=self.pco2)
        self.pH   = self.zerovec + co.pH
        self.DIC  = self.zerovec + co.TC
        self.OmAr = self.zerovec + co.OmAr

    def load():
        return stload.load('tpz_nwnd_ncpm.npz')

    @property
    def fvent(self):
        return self.o2pv / self.mldp

    @property
    def dens(self):
        return sw.dens0(self.salt, self.temp)

    @property
    def o2pv(self):
        return wind2pv(self.temp, self.nwnd)

    
    def calc_CO2_flux(self, tpos=None, w14=True):
        """Calculate CO2 air-sea flux"""
        Sc = schmidt(self.temp[tpos], "co2")
        dpCO2 = self.pco2[tpos] - self.pco2_atm[tpos]
        k0 = co2sys.constants.calculate_K0(self.salt[tpos], self.temp[tpos])
        if w14:
            kappa = 0.251 * self.nwnd[tpos]**2 * (Sc/660)**-0.5
        else:
            kappa = 0.24 * 0.0280 * self.nwnd[tpos]**3 * (Sc/660)**-0.5
        co2fl = kappa * k0 * dpCO2 # (mmol m^-2 d^-1)
        return co2fl
        
    def calc_wwfl(self, wtlen=60):   
        end = self.pv.shape[0]
        wpv = self.pv.copy()
        weights = np.ones(self.pv[wtlen:,...].shape)
        wghtsum = np.ones(self.pv[wtlen:,...].shape)
        for t in np.arange(1,wtlen):
            weights = (weights * (1-self.fvent[wtlen-t+1:end-t+1,...]))
            wghtsum += weights
            wpv[wtlen:,...] += self.pv[wtlen-t:end-t,...] * weights
        wpv[wtlen:,...] = wpv[wtlen:,...] / wghtsum
        wpv[:wtlen,...] = np.nan
        self.wpv = wpv
        self.wwfl = wpv * (self.o2ar) * self.o2st/1000 * self.dens

    def step_co2(self,tpos):
        """Calculate the carbonate system at tpos+1"""
        co = co2sys.CarbonateSystem(self.salt[tpos], self.temp[tpos],
                                 TA=self.talk[tpos], TC=self.DIC[tpos])
        self.pco2[tpos] = co.pCO2
        #print("%03i, %2.2f, %s" % (tpos, self.temp[tpos], co))

        co2fl = self.calc_CO2_flux(tpos=tpos) * self.dens[tpos] / 1000
        self.DIC[tpos+1] = (self.DIC[tpos]
                                - co2fl / self.mldp[tpos] -
                            self.ncpc[tpos])
        co = co2sys.CarbonateSystem(self.salt[tpos],    self.temp[tpos],
                                 TA=self.talk[tpos], TC=self.DIC[tpos+1])
        self.pH[tpos+1]   = co.pH
        self.pco2[tpos+1] = co.pCO2
        self.OmAr[tpos+1] = co.OmAr
        
    def step_o2(self, tpos):
        # Independent vars
        self.o2st[tpos] = sw.satO2(self.salt[tpos], self.temp[tpos])
        o2fl = self.o2pv[tpos] * (self.o2ct[tpos]-self.o2st[tpos])
        self.o2ct[tpos+1] = (self.o2ct[tpos] -
            (o2fl - self.ncpm[tpos])/self.mldp[tpos])

    """
    def step_o2ar(self, tpos):
        # Independent vars
        self.o2ct[self.o2armsk==0] = self.o2st[self.o2armsk==0] 
        nrm = self.o2st[t] / self.arst[t]
        self.o2fl[t] = self.pv[t] * (self.o2ct[t]-self.o2st[t])
        self.arfl[t] = self.pv[t] * (self.arct[t]-self.arst[t])        
        self.o2ct[t+1] = (self.o2ct[t] -
                        (self.o2fl[t] - self.bubl[0]*20 - self.ncpm[t])/self.mldp[t])
        self.arct[t+1] = self.arct[t] - (self.arfl[t] - self.bubl[t])/self.mldp[t]
        self.o2ar = (self.o2ct/self.o2st)/(self.arct/self.arst)*1-1
    """


    """
    def calc_nrm(self, ncpcut=3):
        if ncpcut >= 0:
            self.ncpm[self.ncpm<=ncpcut] = np.nan
        self.ncpm[self.o2ar<=0] = np.nan
        self.nn10[np.isnan(self.ncpm)] = np.nan
        if ncpcut >= 0:
            self.wwfl[self.ncpm<=ncpcut] = np.nan
        self.wwfl[self.o2ar<=0] = np.nan

        self.nrm10 = nodimb(self.wwfl,self.nn10)
        self.nrm01 = nodimb(self.wwfl,self.ncpm)
        self.nrm10[self.nrm10==-1] = np.nan
        self.nrm01[self.nrm01==-1] = np.nan
        self.nrm10[self.nrm10>10] = np.nan
        self.nrm01[self.nrm01>10] = np.nan
        self.nrm10[np.isinf(self.nrm10)] = np.nan
        self.nrm01[np.isinf(self.nrm01)] = np.nan
    """


    def temp_drop(self):
        self.setup_physics()
        self.setup_o2ar()
        self.setup_carb()
        self.temp[100:] = self.temp[100:] - 3
        self.run()
        

    def run(self, label=None):

        #self.salt[100:] = self.salt[100:] - 2
        if HAS_CLICK:
            label= "Run model" if label is None else label
            with click.progressbar(self.tvec[:-1], label=label) as bar:
                for tpos in bar:
                    self.step_co2(tpos)
                    self.step_o2(tpos)
        else:
            for tpos in self.tvec[:-1]:
                self.step_co2(tpos)
                self.step_o2(tpos)

    def to_csv(self, filename, extra_cols={}):
        if type(extra_cols) is str:
            extra_cols = {extra_cols:getattr(self, extra_cols)}
        elif type(extra_cols) is list:
            extra_cols = {key:getattr(self, key) for key in extra_cols}

        df = pd.DataFrame(
            {**{"salt":self.salt, "temp":self.temp, "wind":self.nwnd,
                "pco2_atm":self.pco2_atm,
                "pco2":self.pco2, "talk":self.talk, "pH":self.pH,
                "DIC":self.DIC, "OmAr":self.OmAr}, **extra_cols})
        df.to_csv(filename)


    def plot_temps(self):
        def run(dt):
            self.reset()
            self.temp[100:] = self.temp[100:] - dt
            self.run()
            pl.plot(self.DIC, label="%i $\degree$C drop" % dt, alpha=0.5)
        pl.clf()
        run(1)
        run(2)
        run(3)
        pl.xlabel("Time (days)", fontsize=12)
        pl.ylabel("DIC (mmol/m$^3$)", fontsize=12)
        pl.tick_params(axis='both', which='major', labelsize=12)
        pl.legend()
        pl.title("Windspeed = 10 m/s")

    def plot_winds(self):
        def run(dt):
            self.reset()
            self.temp[100:] = self.temp[100:] - 3
            self.nwnd[:]    = dt
            self.run()
            pl.plot(self.DIC, label="%i m/s wind" % dt, alpha=0.5)
        pl.clf()
        run(5)
        run(10)
        run(15)
        run(20)
        pl.xlabel("Time (days)", fontsize=12)
        pl.ylabel("DIC (mmol/m$^3$)", fontsize=12)
        pl.tick_params(axis='both', which='major', labelsize=12)
        pl.legend()
        pl.title("Temperature drop = 3 $\degree$C")


    def plot_hurricane(self):
        self.reset()
        self.temp[100:103] = self.temp[100:103] - 3
        self.nwnd[100] = 35
        self.run()
        pl.clf()
        #pl.plot(self.DIC, alpha=0.5)
        #pl.ylabel("DIC (mmol/m$^3$)", fontsize=12)
        #pl.ylim(1900,1925)
        #
        #pl.plot(self.pco2, alpha=0.5)
        #pl.ylabel("pCO$_2$ ($\mu$Atm)", fontsize=12)
        #
        pl.plot(self.OmAr, alpha=0.5)
        pl.ylabel("$\Omega_{Ar}$", fontsize=12)
        
        pl.xlabel("Time (days)", fontsize=12)
        pl.tick_params(axis='both', which='major', labelsize=12)
        pl.title("Wind 35 m/s, Temperature drop = 3 $\degree$C")


        
    def plot(self, var1, var2):

        ax1 = pl.gca()
        ax2 = ax1.twinx()

        p1, = ax2.plot(self.DIC,  'r-', label="DIC")
        p2, = ax2.plot(self.pco2, 'g-', label="pCO$_2$")
        ax1.set_xlabel("Time (days)")
        ax1.set_ylabel("DIC")
        ax2.set_ylabel("pCO$_2$")
        """
        lines = [p1, p2, p3]
        host.legend(lines, [l.get_label() for l in lines])

        for ax in [par1, par2]:
            ax.set_frame_on(True)
            ax.patch.set_visible(False)

        plt.setp(ax.spines.values(), visible=False)
        ax.spines["right"].set_visible(True)
        """
        ax1.yaxis.label.set_color(p1.get_color())
        ax2.yaxis.label.set_color(p2.get_color())
        #par1.spines["right"].set_edgecolor(p2.get_color())
        #par2.spines["right"].set_edgecolor(p3.get_color())
        ax1.tick_params(axis='y', colors=p1.get_color())
        ax2.tick_params(axis='y', colors=p2.get_color())

        

#####################################

    
def run_model(st):
    for t in np.arange(0,364):
        # Independent vars
        self.o2ct[self.o2armsk==0] = self.o2st[self.o2armsk==0] 
        nrm = self.o2st[t]/self.arst[t]
        self.o2fl[t] = self.pv[t] * (self.o2ct[t]-self.o2st[t])
        self.arfl[t] = self.pv[t] * (self.arct[t]-self.arst[t])        
        self.o2ct[t+1] = (self.o2ct[t] -
            (self.o2fl[t] - self.bubl[0]*20 - self.ncpm[t])/self.mldp[t])
        self.arct[t+1] = (self.arct[t] -
            (self.arfl[t] - self.bubl[t])/self.mldp[t])
    self.o2ar = (self.o2ct/self.o2st)/(self.arct/self.arst)*1-1
    return st


#####################################

    
def steady():    
    st = setup()
    st = setup_o2ar(st)
    st = run_model(st)
    plot1D(st)

def pulse():
    st = setup()
    self.nwnd[100:110] = 15
    self.ncpm[100:110] = 50
    st = setup_o2ar(st)
    st = run_model(st)
    pl.figure(1)
    plot1D(st)
    pl.figure(2)
    plotReuer(st)
    return st

def growth_season(tpp,im,tl=180,i1=None,i2=None,j1=None,w=False):
    if w:
        tp = tpp
        st = setup(tp.nwnd[:,10:12,10:12].shape)
        self.nwnd[...] = w
        self.ncpm = np.zeros(self.nwnd.shape) *0.01
    else:
        tp = subset(tpp,im)
        st = setup(tp.nwnd[:,:im,:].shape)
        self.nwnd = tp.nwnd[:,...].copy()
        self.ncpm = np.zeros(self.nwnd.shape) *0.01
    xvec = np.arange(tl)
    tsin = np.zeros(self.nwnd.shape[0]) *0.01
    tsin[90:90+tl,...] = np.sin(xvec/((len(xvec)-1)/np.pi))*40
    self.ncpm[:,...] = tsin[:,np.newaxis,np.newaxis]
    self.ncpf = self.ncpm.copy()

    st = setup_o2ar(st)
    st = run_model(st)
    st = calc_wwfl(st)
    st = calc_nn10(st)
    st = calc_nrm(st,ncpcut=-1)
    return st

def legtext(txt,nrm):
    return ('%s  avg=%.2f   std=%.2f' %
            (txt,np.mean(nrm[~np.isnan(nrm)]),
             np.std(nrm[~np.isnan(nrm)])) )

def season_hist(tp,im):
    pl.close(1)
    F = pl.figure(1,(11,8))
    slvec =  [30,90,150,210]
    pl.clf()
    legvec = []
    for n,s in enumerate(slvec):
        st = growth_season(tp,im=im,tl=s)
        pl.hist(self.nrm10[~np.isnan(self.nrm10)],np.arange(-2,2,0.02),
                histtype='step')
        legvec.append(legtext('%i days' % s,self.nrm10))
    pl.legend(legvec)
    pl.title('Misfit for different lenghts of Growth Seasons')
    pl.ylabel('Number of data points')
    pl.xlabel('Misfit ( (bioflux-NCP)/ncp )')
    pl.savefig('figs/toymodel/seasonhiself.png')

def season_wind(tp,im,i1=10,j1=10):
    slvec =  [30,90,150,210]
    pl.close(1)
    F = pl.figure(1,(11,8))
    w = 10
    for n,s in enumerate(slvec):
        print(n)
        st = growth_season(tp,im,tl=s,w=w)
        pl.subplot(4,1,n+1)
        pl.plot(self.ncpf[:,0,0])
        pl.plot(self.nn10[:,0,0])
        pl.plot(self.wwfl[:,0,0])
        pl.xticks(np.arange(0,400,50),[])
        pl.yticks([10,30])
        pl.title('Growth period = %03i days' % s)
        pl.legend ( ('NCP instant','NCP 10d avg', 'Bioflux'))
        pl.xlim(0,365)
        pl.ylabel('mmol/m2')
    pl.xticks(np.arange(0,400,50),np.arange(0,400,50))
    pl.xlabel('Days')
    pl.suptitle('Constant wind (%i m/s)' % w)
    pl.savefig('figs/toymodel/seasonwind_%03i.png' % w)
    
    slvec = np.arange(30,250,5)
    wvec = [5,10,15]

    pl.close(2)
    F = pl.figure(2,(11,8))
    legvec = []
    for w in wvec:
        msvec = []
        for s in slvec:
            st = growth_season(tp,im,tl=s,w=w)
            msvec.append(np.mean(abs(self.nrm10[~np.isnan(self.nrm10)])))
        pl.plot(slvec,msvec,'.')
        legvec.append('%02i m/s' % w)
    pl.legend(legvec)
    pl.ylim(0,1)
    pl.ylabel('Average misfit during growth season, sign removed.')
    pl.xlabel('Growing season (days)')
    pl.suptitle('Constant wind speed vs length of growth season')
    pl.savefig('figs/toymodel/seasonwind_all.png')
  

def season_plot(tp,i1=10,j1=10):
    slvec =  [30,90,150,210]
    pl.figure(2)
    pl.clf()
    for n,s in enumerate(slvec):
        print(n)
        st = growth_season(tp,tl=s,i1=i1,j1=j1)
        pl.figure(2)
        pl.subplot(4,1,n+1)
        pl.plot(self.nn10)
        pl.plot(self.wwfl)
        pl.xticks(np.arange(0,400,50),[])
        pl.yticks([10,30])
        pl.title('Growth period = %03i days' % s)
    pl.xticks(np.arange(0,400,50),np.arange(0,400,50))
    pl.xlim(0,365)
    
    slvec = np.arange(30,250,5)
    msvec = []
    for s in slvec:
        st = growth_season(tp,tl=s,i1=i1,j1=j1)
        msvec.append(np.mean(abs(self.msft)))


    pl.legend ( ('NCP 10d avg', 'Bioflux'))
    pl.figure(3)
    #pl.clf()
    pl.plot(slvec,msvec,'.')
    #pl.xlim(20,200)
    pl.ylim(0,1)

    pl.ylabel('Average misfit during growing season')
    pl.xlabel('Growing sweason (days)')
    pl.figure(2)
    pl.savefig('figs/seasonplot_%03i_%03i_2.png' % (i1,j1) )
    pl.figure(3)
    pl.savefig('figs/seasonplot_%03i_%03i_3.png' % (i1,j1) )

def season_plot_2D(tpp,im=30):
    class tp: pass
    tpvars = ['llat', 'llon', 'mldp','ncpm','nwnd', 'temp']    
    for v in tpvars:
        tp.__dict__[v] = tpp.__dict__[v][...,:im,:].copy()
    slvec =  [30,90,150,210]
    pl.figure(2)
    pl.clf()
    for n,s in enumerate(slvec):
        st = growth_season(tp,tl=s)
        pl.figure(n+1)
        ncphist(st,self.msft)

def realwind(tp):
    st = setup(im,jm)

    self.nwnd = tp.nwnd[:,40,10]
    #self.ncpm = tp.ncpm[:,40,10]
    #self.temp = tp.temp[:,40,10]

    self.o2ct[0] = 200
    self.arct[0] = 11
    self.o2st = sw.satO2(self.salt, self.temp)
    self.arst = sw.satAr(self.salt, self.temp)
    self.pv   = wind2pv(self.temp,self.nwnd)
    
    st = run_model(st)
    pl.figure(1)
    plot1D(st)
    pl.figure(2)
    plotReuer(st)
    return st

def realwind2D(tpp,im,realvars=['nwnd','ncpm']):
    tp = subset(tpp,im)
    st = setup(tp.nwnd[:,:im,:].shape)
    if 'nwnd' in realvars:
        self.nwnd = tp.nwnd[:,:im,:].copy()
    if 'ncpm' in realvars:
        self.ncpm = tp.ncpm[:,:im,:].copy()
        self.ncpm[np.isnan(self.ncpm)] = 0
    if 'mldp' in realvars:
        self.mldp = tp.mldp[:,:im,:].copy()
    if 'temp' in realvars:
        self.temp = tp.temp[:,:im,:].copy()
    if 'o2armsk' in realvars:
        self.o2armsk = tp.o2ar[:,:im,:].copy()
        self.o2armsk[self.o2armsk>0] = 1
        self.o2armsk[self.o2armsk<0] = 0
    st = setup_o2ar(st)
    st = run_model(st)
    st = calc_wwfl(st)
    st = calc_nn10(st)
    st = calc_nrm(st)
    return st

def hist_realwind2D(tp,im,realvars=['nwnd','ncpm']):
    st = realwind2D(tp,im,realvars=realvars)

    pl.clf()
    nvec = self.nrm10[~np.isnan(self.nrm10)]
    h = pl.hist(nvec,np.arange(-2,2,0.02),histtype='step')
    nvec = self.nrm01[~np.isnan(self.nrm01)]
    h = pl.hist(nvec,np.arange(-2,2,0.02),histtype='step')

    if 1==0:
        tp = calc_nrm(tp)
        self.nrm01[np.isnan(tp.nrm01)] = np.nan
        nvec = tp.nrm10[~np.isnan(tp.nrm10)]
        h = pl.hist(nvec,np.arange(-2,2,0.02),histtype='step')
        nvec = tp.nrm01[~np.isnan(tp.nrm01)]
        h = pl.hist(nvec,np.arange(-2,2,0.02),histtype='step')
    pl.legend( ('NCP 10-day mean','NCP instant','full10d','fullinst') )
    return st


def ncphist(st):

    self.nrm10[np.isinf(self.nrm10)]=np.nan
    nvec= [0,10,20,30,40,50]
    leglist = []
    def hist(nrm):
        return pl.hist(nrm[~np.isnan(nrm)],np.arange(-2,2,0.02),
                       histtype='step')
    hist(self.nrm10)
    legliself.append(legtext('All data',self.nrm10))
    for n in nvec:
        nrm = self.nrm10.copy()
        nrm[self.nn10>n+10] = np.nan
        nrm[self.nn10<n]    = np.nan 
        hist(nrm)
        legliself.append(legtext('%02i-%02i' % (n,n+10),nrm))
    pl.legend (leglist)
    pl.title('Misfit for different NCP values')
    pl.ylabel('Number of data points')
    pl.xlabel('Misfit ( (bioflux-NCP)/ncp )')
    pl.savefig('figs/toymodel/ncphiself.png')


def plot1D(st):
    pl.clf()
    pl.subplot(4,1,1)
    pl.plot(self.o2ct/self.o2st*100 - 100)
    pl.plot(self.arct/self.arst*100 - 100)
    pl.ylabel('Percent')
    lg = pl.legend(('O2','Ar'))
    lg.get_frame().set_fill(False)
    pl.title('Saturation')
    pl.ylim(-10,10)
    pl.xticks([])

    pl.subplot(4,1,2)
    pl.plot(self.o2fl)
    pl.plot(self.arfl*20)
    lg = pl.legend(('O2','Ar'))
    lg.get_frame().set_fill(False)
    pl.title('Flux')
    pl.ylabel('mmol / d')
    pl.ylim(0,100)
    pl.xticks([])
    
    pl.subplot(4,1,3)
    pl.plot(self.o2ar*100)
    pl.title('O2Ar Saturation')
    pl.ylim(-10,10)
    pl.ylabel('Percent')
    pl.xticks([])
    
    pl.subplot(4,1,4)
    pl.plot(self.o2ar * self.pv*self.o2st)
    pl.plot(self.ncpm)
    pl.title('O2bio Flux')
    pl.ylabel('mmol / d')
    pl.ylim(0,100)

def plotReuer(st, wtlen=60):
    st = calc_wwfl(st)
    pl.clf()
    pl.subplot(4,1,1)
    pl.plot(self.nwnd)
    pl.plot(self.pv)
    lg = pl.legend(('Wind','PV'))
    lg.get_frame().set_fill(False)
    pl.title('Wind (m/s) & PV (m/d)')
    pl.ylim(0,20)
    pl.xticks([])
    pl.xlim(0,400)
    
    pl.subplot(4,1,2)
    pl.plot(self.fvent)
    pl.title('F-vent')
    pl.ylabel('m^2 / d')
    pl.ylim(0,1)
    pl.xlim(0,400)
    
    pl.subplot(4,1,3)
    pl.plot(self.wpv)
    pl.title('Weighted PV')
    pl.ylabel('m / d')
    pl.ylim(0,20)
    #pl.ylabel('Percent')
    pl.xlim(0,400)
    
    pl.subplot(4,1,4)
    pl.plot(self.wwfl)
    pl.plot(self.ncpm)

    pl.title('O2bio Flux')
    pl.ylabel('mmol / d')
    pl.ylim(0,100)
    pl.xlim(0,400)


#md = BoxModel()
#md.run()             
