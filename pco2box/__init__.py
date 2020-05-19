"""Model to study the solubility and gasexchange in a generalized ocean box """
import os

from .timeseries import plot_perturbation
from .recovery import plot_tmat_wind_mld


FIGDIR = os.path.abspath("./pCO2box_figs")
if not os.path.exists(FIGDIR):
    os.makedirs(FIGDIR)
