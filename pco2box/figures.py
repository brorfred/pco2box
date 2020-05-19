
import numpy as np
import pylab as pl
import seaborn as sns
#import cartopy.crs as ccrs
#crs = ccrs.PlateCarree()
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter

jdvec = np.arange(735599, 736481)[:730]
months = MonthLocator(range(1, 13), bymonthday=1, interval=3)
monthsFmt = DateFormatter("%b")

#def east_map():
#    pl.clf()
#    fig = pl.figure()
#    ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal(
#        central_longitude=-70, central_latitude=45))
#    ax.set_extent([-81, -58, 23, 47], crs=ccrs.PlateCarree())
#    ax.coastlines('10m',zorder=3, color="0.5")
#
#    latvec = np.load("indata/merc_uscoast_latlon.npz")["lat"][[300,400,460]]
#    lonvec = np.load("indata/merc_uscoast_latlon.npz")["lon"][[710,770,840]]
#    pl.scatter(lonvec, latvec, transform=crs)
#    transform = crs._as_mpl_transform(pl.gca())
#    for txt,lon,lat in zip(["SAB","MAB", "NWA"], lonvec, latvec):
#        pl.annotate(txt, xy=(lon,lat), xycoords=transform)
#    pl.savefig("boxmodel_eastwest_east_map.png", bbox_inches="tight")

def fldplot(md, fldname, ax=None, **kwargs):
    """Plot list from model"""
    if ax is None:
        ax = pl.gca()
    for args in {"lw":1, "alpha":0.5}.items():
        kwargs[args[0]] = kwargs.get(args[0], args[1])
        
    p = ax.plot(jdvec,getattr(md, fldname), **kwargs)
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(monthsFmt) 
    ax.xaxis.set_major_locator(months) 
    ax.grid(True)
    ax.set_xlim("2015-01-01","2017-01-01")
    return p[0].get_color()

def plot_timeseries(md, **kwargs):
    tsp=3
    sns.set_style("whitegrid")
    clrs =sns.color_palette()
    if not "axes" in kwargs:
        fig,axes = pl.subplots(tsp,1, sharex=True, num=1,figsize=(6,6))
    else:
        axes = kwargs.pop("axes")
    pl.subplots_adjust(hspace=0.03)

    def line(ax, fldname, ylabel, ylim=(None,None),axpos="left", yticks=None,
                 **kwargs):
        if axpos == "right":
            kwargs["c"] = clrs[1]
            ax = pl.twinx(ax=ax)
            ax.yaxis.label.set_color(clrs[1])
            ax.tick_params(axis='y', colors=clrs[1])
        else:
            kwargs["c"] = clrs[0]
            ax.yaxis.label.set_color(clrs[0])
            ax.tick_params(axis='y', colors=clrs[0])
        fldplot(md, fldname, ax=ax, **kwargs)
        if fldname =="DIC":
            kwargs["label"] = None
            fldplot(md, "DICeq", ax, ls=":", **kwargs)
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        if yticks is not None:
            ax.set_yticks(yticks)

    args = [axes[0], "temp", "Temperature $\degree$C", [0,30]]
    line(*args, yticks=[7.5,15,22.5], **kwargs)
    args = [axes[0], "salt", "Salinity ", [32,36]]
    line(*args, axpos="right", yticks=[33,34,35], **kwargs)
    args = [axes[1], "pco2", "pCO$_2$ $\mu{}atm$", [300,1100]]
    line(*args, yticks=[500,700,900], **kwargs)
    args = [axes[1], "pH", "pH", [7.5,8.5], ]
    line(*args, axpos="right", yticks=[7.75, 8.0, 8.25], **kwargs)
    args = [axes[2], "DIC",  "DIC mmol m$^{-3}$", [1900,2300]]
    line(*args, yticks=[2000,2100,2200], **kwargs) 
    args = [axes[2], "OmAr",  "$\Omega$", [0.6, 3.0]]   
    line(*args, axpos="right", yticks=[1.2, 1.8, 2.4],**kwargs)
