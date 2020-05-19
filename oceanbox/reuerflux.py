import numpy as np

import seawater as sw
from . import wind2pv


def pv2wpv(temp, mld, wind, wtlen=60, salt=False, dens=None):
    """Calculate bio-o2flux using weighted wind """
    dens = sw.dens0(salt, temp) if dens is None else dens
    pv    = wind2pv(temp,wind)
    fvent = pv / mld
    fvent[np.isnan(fvent)] = 0
    pv[np.isnan(pv)] = 0
    if wtlen == 0:
        wpv = pv
    else:
        end=pv.shape[0]
        wpv = pv.copy()
        weights = np.ones(pv[wtlen:,...].shape)
        wghtsum = np.ones(pv[wtlen:,...].shape)
        for t in np.arange(1,wtlen):
            weights = (weights * (1-fvent[wtlen-t+1:end-t+1,...]))
            wghtsum += weights
            wpv[wtlen:,...] += pv[wtlen-t:end-t,...] * weights
        wpv[wtlen:,...] = wpv[wtlen:,...] / wghtsum
    wpv[:wtlen,...] = np.nan
    return wpv

def wpvweight(temp, mld, wind, wtlen=60, salt=False, dens=False):
    """Calculate bio-o2flux using weighted wind """
    if type(dens) == bool : dens = sw.dens0(salt, temp)
    pv    = wind2pv(temp,wind)
    fvent = pv / mld
    fvent[np.isnan(fvent)] = 0
    pv[np.isnan(pv)] = 0
    if wtlen == 0:
        wpv = pv
    else:
        end=pv.shape[0]
        wpv = field.copy()
        weights = np.ones(pv[wtlen:,...].shape)
        wghtsum = np.ones(pv[wtlen:,...].shape)
        for t in np.arange(1,wtlen):
            weights = (weights * (1-fvent[wtlen-t+1:end-t+1,...]))
            wghtsum += weights
            wpv[wtlen:,...] += pv[wtlen-t:end-t,...] * weights
        wpv[wtlen:,...] = wpv[wtlen:,...] / wghtsum
    wpv[:wtlen,...] = np.nan
    return wpv[:,0,...]

def test_range(patm):
    if patm.min() < 500 or patm.min() > 1500: 
        print("Atmospheric pressure on in range")
        raise

def pv2fl(o2ar,temp,mld,wind,wtlen=60,
          salt=False, o2st=False, dens=False, patm=False):
    """Calculate bio-o2flux using weighted wind """
    if type(o2st) == bool : o2st = sw.satO2(salt, temp)
    if type(patm) == np.ndarray:
        print("Using atm press in reuerflux")
        test_range(patm)
        o2st=o2st * (patm/1013.25)
    if type(dens) == bool : dens = sw.dens0(salt, temp)
    pv    = wind2pv(temp,wind)
    fvent = pv / mld
    fvent[np.isnan(fvent)] = 0
    pv[np.isnan(pv)] = 0
    if wtlen == 0:
        wpv = pv
    else:
        end=pv.shape[0]
        wpv = pv.copy()
        weights = np.ones(pv[wtlen:,...].shape)
        wghtsum = np.ones(pv[wtlen:,...].shape)
        for t in np.arange(1,wtlen):
            weights = (weights * (1-fvent[wtlen-t+1:end-t+1,...]))
            wghtsum += weights
            wpv[wtlen:,...] += pv[wtlen-t:end-t,...] * weights
        wpv[wtlen:,...] = wpv[wtlen:,...] / wghtsum
    wpv[:wtlen,...] = np.nan
    return (wpv[:,np.newaxis] * (o2ar/100) * o2st/1000 * dens)[:,0,...]


def bottle2o2ar(o2ar_air, temp, salt):

    air = 20.946/0.9340
    o2ar_sat = (sw.satO2(salt, temp) / sw.satAr(salt, temp)/air-1) * 1000
    return o2ar_air / o2ar_sat - 1
