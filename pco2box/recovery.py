import os

import numpy as np
import pylab as pl
 
from .model import perturbation
from .figures import plot_timeseries

FIGDIR = os.path.abspath("./pCO2box_figs")


def recovery_time(wind=7, mldp=30, temp=12, salt=33, pfunc=perturbation):
    """Calculate recovery time for a specific mld and windspeed
    
    Run the model with and without a pertubation calculate the time it
    takes for the difference to be less than 1/e of the initial difference.

    Parameters
    ----------
    wind : float, optional
        Prescribed fixed wind
    mldp : float, optional
        Prescribed fixed MLD
    temp : float, optional
        Prescribed fixed temperature
    salt : float, optional
        Prescribed fixed salinity
    pfunc : function, optional
        Function to calculate pertubation and baseline

    Returns
    -------
    time : int
        Time in days to return to 1/e difference between pertubation
        and baseline timeseries.
    """
    mds,mdp = pfunc(wind=wind, mldp=mldp, temp=temp, salt=salt)
    pvec = (mdp.DIC-mds.DIC)[mpos+1:]/(mdp.DIC-mds.DIC)[mpos+1] 
    pvec[pvec<1/np.exp(1)] = 0 
    return np.nonzero(pvec==0)[0].min()

def perturbation_recovery_matrix(parname="DIC"):
    """Compute array with recovery timescales for different winds and mld's
    
    Returns an array with recovery timescales based on different 
    Mixed Layer Depths (MLDs) and windspeeds defined by the matrices

        mmat,wmat = np.mgrid[10:100:5,1:20:2]

    The boxmodel is run with fixed a prescribed temperature of temp 12°C and a
    prescribed salinity of 33 PSU. The pertubation is set to occur after
    165 day. This values has no effect on the results since T and S are fixed.

    use the function 'plot_tmat_wind_mld' to visualize the results.

    Parameters
    ----------
    parname : string, optional
        Name of carbon species to base the recovery (default is DIC)

    Returns
    -------
    recovery_time_matrix : ndarray
        matrix with wind as xaxis and MLD as y-axis
    """
    mpos = 165 # Time of pertubation. Has no effect since T and S are fixed
    mmat,wmat = np.mgrid[10:100:5,1:20:2]
    tmat = mmat * np.nan
    for m,mldp in enumerate(mmat[:,0]):
        for n,wind in enumerate(wmat[0,:]):
            mds,mdp = perturbation(wind=wind, mldp=mldp, mpos=mpos)
            tvecs = getattr(mds, parname)
            tvecp = getattr(mdp, parname)
            pvec = (tvecp-tvecs)[mpos+1:]/(tvecp-tvecs)[mpos+1]             
            pvec[pvec<1/np.exp(1)] = 0 
            try:
              tmat[m,n] = np.nonzero(pvec==0)[0].min()
            except:
                print(mldp,wind)
                pass
    return tmat

def plot_tmat_wind_mld(tmat=None, parname="DIC", color="k", text=""):
    """Plot tmat array"""
    mmat,wmat = np.mgrid[10:100:5,1:20:2]
    if tmat is None:
        tmat = perturbation_recovery_matrix(parname=parname)
    CS = pl.contour(wmat, mmat, tmat,
                    [10,20,30,40,50,60,80,100,150,200],
                    colors=color)
    pl.gca().clabel(CS, inline=1, fontsize=10, fmt="%i")
    CS.collections[0].set_label(text)
    pl.xlim(3,18)
    pl.xlabel("Wind speed (m/s)")
    pl.ylabel("Mixed Layer Depth (m)")
    pl.title("Time to recover from 51 mmol DIC perturbation (days)")
    pl.savefig(os.path.join(FIGDIR, f"recovery_matrix.pdf"))

    return CS

#In [87]: plot([0.01, 0.3, 0.4, 0.5, 1, 1.5,2, 3, 4, 5], [30,28,27,26,22,18,14, 12, 10, 9], ".-")                      
