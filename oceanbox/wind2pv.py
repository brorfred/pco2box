import warnings
import numpy as np

def wind2pv(temp, wind, model='swee', gas="o2"):
    temp = np.array(temp)
    wind = np.array(wind)
    if temp.ndim > wind.ndim:
        temp = temp[:,0,...]
    assert not wind[~np.isnan(wind)].min() < 0, "Negative wind speeds."
    if any([s in model.lower() for s in ['wann','w92']]):
        return _wind2pv_wann1992(temp, wind, gas)
    elif  any([s in model.lower() for s in ['swee','sw07']]):
        return _wind2pv_sweeney2007(temp, wind, gas)

def _wind2pv_wann1992(temp, wind, gas="o2"):
    sc = schmidt(temp, gas)
    return (0.31 * wind**2 * (sc/660)**(-0.5)) * (24./100) 

def _wind2pv_sweeney2007(temp, wind, gas="o2"):
    sc = schmidt(temp, gas)
    return (0.27 * wind**2 * (sc/660)**(-0.5)) * (24./100) 

#def _wind2pv_


def schmidt(temp, gas="o2"):
    """Calculate the Schmidt number for different gases
    
    Formulation from Saramiento, Chemical Ocenaography, p 85. 
    Based on Wilke and Chang (1955) via Wanninkhof (1992).
    """
    if gas =="o2":
        Sc = {"A":1953.4, "B":128.00, "C":3.9918, "D":0.050091}
    elif gas == "ar":
        Sc = {"A":1909.1, "B":125.09, "C":3.9012, "D":0.048953}
    elif gas == "co2":
        Sc = {"A":2073.1, "B":125.62, "C":3.6276, "D":0.043219}
    return Sc["A"] - Sc["B"]*temp + Sc["C"]*temp**2 - Sc["D"]*temp**3

