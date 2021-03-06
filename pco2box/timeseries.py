import os

import numpy as np
import pylab as pl
 
from .model import run_model, perturbation
from .figures import plot_timeseries

FIGDIR = os.path.abspath("./pCO2box_figs")

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
    pl.savefig(os.path.join(FIGDIR, f"pertubation_timeseries_{pref}.pdf"))

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
    
