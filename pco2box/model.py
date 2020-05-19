
import os

import numpy as np
import pylab as pl
import pandas as pd
from scipy.stats import linregress
 
import oceanbox
import co2sys

from .figures import plot_timeseries

DATADIR = os.path.join(os.path.dirname(__file__), "indata")

def calc_talk(salt, pref):
    """Calculate Total Alkalinity
    
    Return Total Alkalinity based on parametrizations for different parts
    of the US East Coast. Based on personal communications with Wei-Jun Cai
    and [REF]
    
    Parameters
    ----------
    salt : array_like
        Salinity used to calculate TALk
    pref : string
        Part of the coast [se=South Atlantic Bight, 
                           me=Mid Atlantic Bight, 
                           ne=North Atlantic outside Gulf of Maine]

    Returns
    -------
        TAlk : array_like
            Calculated Total Alkalinity
    """
    if pref == "se":
        talk = 48.7 * salt + 608.8
    if pref == "me":
        talk = 46.6 * salt + 670.6
        #md.salt = md.salt - 1
    if pref == "ne":
        talk = 39.1 * salt + 932.7
    return talk

def setup_model(pref="ne", wind=7, mldp=30, wintermixing=True, winterwind=True):
    """Setup Box model
    
    Initiate the boxmodel class with data from the indata/mercator_tseries.h5 
    file. Wind and Mixed Layer Depths (MLDs) can be prescribed to fixed values 
    and winter conditions can be turned on or off.

    Parameters
    ----------
    pref : string, optional
        Part of the coast [se=South Atlantic Bight, 
                           me=Mid Atlantic Bight, 
                           ne=North Atlantic outside Gulf of Maine]
                           ca=California Coastal Current]
    wind : float, optional
        Prescribed fixed wind
    mldp : float, optional
        Prescribed fixed MLD
    wintermixing : bool, optional
        Simulate wintermixing by removing 0.1 mM DIC daily during winter
    winterwind : bool, optional
        Simulate winter condtions by increasing winds 1.5 times during winter
    """
    df = pd.read_hdf(os.path.join(DATADIR, "mercator_tseries.h5"))
    md = oceanbox.BoxModel(svec=730)
    md.nwnd = np.zeros((365)) + wind
    if winterwind:
        md.nwnd[:60]  = wind*1.5
        md.nwnd[270:] = wind*1.5
    md.nwnd = np.hstack((md.nwnd, md.nwnd))  
    md.mldp[:] = mldp

    def loop_years(vec):
        slope = np.linspace(0,(vec[0] - vec[364]), 365)
        vec = vec + slope
        return np.hstack((vec, vec))

    md.temp = loop_years(df[pref + "temp"][:365])
    md.salt = loop_years(df[pref + "salt"][:365])
    try:
        md.ncpc = loop_years(df[pref + "ncpm"][:365])
    except KeyError:
        md.ncpc = loop_years(np.zeros(365))
    if wintermixing:
        md.ncpc[270:] = md.ncpc[270:] - 0.1
        md.ncpc[:60]  = md.ncpc[:60]  - 0.1
    return md

def run_model(pref, reg, temp=None, salt=None, wind=7, mldp=30,
        deepmld1={"pert":0, "tpos":165},
        deepmld2={"pert":0, "tpos":553}, 
        uppwell1={"pert":0, "tpos":553, "pdays":10},
        uppwell2={"pert":0, "tpos":553, "pdays":10}):
    """Run the boxmodel in pertubation mode

    Setup and run the box model with the ability to apply an pertubation.

    pref : string
        Part of the coast [se=South Atlantic Bight, 
                           me=Mid Atlantic Bight, 
                           ne=North Atlantic outside Gulf of Maine]
                           ca=California Coastal Current]
    reg : string
        Region name to be used in title and file name
    wind : float, optional
        Prescribed fixed wind
    mldp : float, optional
        Prescribed fixed MLD
    temp : float, optional
        Prescribed fixed temperature
    salt : float, optional
        Prescribed fixed salinity
    deepmld1 : dict, optional
        Dict describing a mixed layer deepening pertubation
    deepmld2 : dict, optional
        Dict describing a mixed layer deepening pertubation
    uppwell1 : dict, optional
        Dict describing a upwelling pertubation
    uppwell2 : dict, optional
        Dict describing a upwelling pertubation
    """
    md = setup_model(pref=pref, mldp=mldp, wind=wind)
    if temp is not None:
        md.temp[:] = temp
    if salt is not None:
        md.salt[:] = salt

    def mix_deepmld(deepmld):
        md.ncpc[deepmld["tpos"]] = md.ncpc[deepmld["tpos"]] - deepmld["pert"]
        md.mldp[deepmld["tpos"]:deepmld["tpos"]+10] = mldp*2
        ncpslice = slice(deepmld["tpos"]+5, deepmld["tpos"]+15)

    def mix_upwell(uppwell):
        tpos = uppwell["tpos"]
        pdays = uppwell["pdays"]
        md.ncpc[tpos:tpos+pdays+1] -= uppwell["pert"]/pdays

    if deepmld1["pert"] != 0:
        mix_deepmld(deepmld1)
    if deepmld2["pert"] != 0:
        mix_deepmld(deepmld2)
    if uppwell1["pert"] != 0:
        mix_upwell(uppwell1)
    if uppwell2["pert"] != 0:
        mix_upwell(uppwell2)

    talk = calc_talk(md.salt, "ne")
    co = co2sys.CarbonateSystem(md.salt, md.temp, TA=talk, pCO2=talk*0+395)
    md.DICeq = co.TC            
    md.setup_carb()
    md.run(label="Run %s" % reg)
    return md

def perturbation(wind=7, mldp=30, temp=12, salt=33, mpos=165, pdays=40):
    """Run the model twice, once with perturbation once without
    
    Parameters
    ----------
    pref : string, optional
        Part of the coast [se=South Atlantic Bight, 
                           me=Mid Atlantic Bight, 
                           ne=North Atlantic outside Gulf of Maine]
                           ca=California Coastal Current]
    wind : float, optional
        Prescribed fixed wind
    mldp : float, optional
        Prescribed fixed MLD
    temp : float, optional
        Prescribed fixed temperature
    salt : float, optional
        Prescribed fixed salinity
    mpos : int, optional
        Day of year of pertubation
    pdays : int, optional
        Length of upwelling pertubation

    Returns
    -------
    mds : model class instance
        Background case
    mpd : model class instance
        Pertubation case
    """
    kwargs = dict(pref="ne", reg="NWA",
                  temp=temp, salt=salt, wind=wind, mldp=mldp,
                  deepmld1={"pert":0, "tpos":mpos})
    mds = run_model(**kwargs)
    kwargs["deepmld1"] = {"pert":17*3, "tpos":mpos}
    mdp = run_model(**kwargs)
    return mds,mdp

def plot_perturbation(wind=7, mldp=30, pref="ca"):
    """Run the model in pertubation mode and plot the results"""
    pl.clf()
    fig,axes = pl.subplots(3, 1, sharex=True, num=1,figsize=(6,6))
    model_kws = dict(pref=pref, reg="pert",
                temp=None, salt=None, wind=wind, mldp=mldp,
                deepmld1={"pert":0, "tpos":165},
                deepmld2={"pert":0, "tpos":553},
                uppwell1={"pert":0, "tpos":165, "pdays":5},
                uppwell2={"pert":0, "tpos":553, "pdays":5})
    md = run_model(**model_kws)
    plot_timeseries(md, axes=axes, alpha=0.5) 
    if pref == "ca":
        preftxt = "CCS"
        model_kws["uppwell1"]["pert"] = 82.5
        model_kws["uppwell2"]["pert"] = 165
    else:
        preftxt = "NWA"
        model_kws["deepmld1"]["pert"] = 17
        model_kws["deepmld2"]["pert"] = 34
    md = run_model(**model_kws)
    plot_timeseries(md, axes=axes, alpha=1)
    pl.suptitle(
        f"Perturbations, temp and salt for {preftxt}, wind:{wind}m/s, mld:{mldp}m")
    pl.savefig(f"figs/pertubation_timeseries_{pref}.pdf")

def to_csv(wind=7, mldp=30, pref="ca"):
    """Save model results to csv
    
    Parameters
    ----------
    wind : float, optional
        Prescribed fixed wind
    mldp : float, optional
        Prescribed fixed MLD
    pref : string, optional
        Part of the coast [ne=North Atlantic outside Gulf of Maine]
                           ca=California Coastal Current]    
    """
    ptext = "CCS" if pref == "ca" else "NWA"
    if pref == "ca":
        preftxt = "CCS"
        model_kws = dict(pref=pref, reg="pert",
                     temp=None, salt=None, wind=wind, mldp=mldp,
                     deepmld1={"pert":0, "tpos":165},
                     deepmld2={"pert":0, "tpos":553},
                     uppwell1={"pert":82.5, "tpos":165, "pdays":5},
                     uppwell2={"pert":165,  "tpos":553, "pdays":5})
    else:
        preftxt = "NWA"
        model_kws = dict(pref=pref, reg="pert",
                     temp=None, salt=None, wind=wind, mldp=mldp,
                     deepmld1={"pert":17, "tpos":165},
                     deepmld2={"pert":34, "tpos":553},
                     uppwell1={"pert":0,  "tpos":165, "pdays":5},
                     uppwell2={"pert":0,  "tpos":553, "pdays":5})

    md = run_model(**model_kws)
    md.to_csv(f"boxmodel_baseline_{ptext}__wind_{wind:02}__mldp_{mldp:03}.csv")
    for key in ["deepmld1", "deepmld1", "uppwell1", "uppwell2"]:
        model_kws[key] = {"pert":0, "tpos":165}
    md = run_model(**model_kws)
    md.to_csv(f"boxmodel_pertubation_{ptext}__wind_{wind:02}__mldp_{mldp:03}.csv")

 
