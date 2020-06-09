
# pCO2 box
> Box model simulating carbon solubility and air-sea exchange in the ocean.


This package includes all code neccessary to run simulation and generate figues for the paper "Controls on Surface Water Carbonate Chemistry along North American Ocean Margins" by Cai et al.

## Background and boundary conditions
To study the relationships that control surface ocean carbonate parameters and air-sea CO2 fluxes, we use a box model to representing an idealized mixed layer. The box model is simulating the carbonate systems as a homogenous surface-ocean water mass located either in the North-West Atlantic (NWA, 40.1-42.5°N, 69.2-67.5°W), South Atlantic Bight (28.3°N–29.7°N, 77.6°W–79.2°W), or the California Current System (CCS, 40.8–41.6°N, 124.5–125°W). The model is driven by realistic temperatures and salinities from the Mercator 1/12° data-assimilated General Circulation Model with a daily resolution in time[^1]. 

Biological production is assessed from weekly averages of daily changes in satellite derived Chlorophyll for the year 2015 using the daily MODIS aqua 4-km Level 3 product (NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group. Moderate-resolution Imaging Spectroradiometer (MODIS) Aqua Chlorophyll Data; 2018 Reprocessing. NASA OB.DAAC, Greenbelt, MD, USA. doi: data/10.5067/AQUA/MODIS/L3M/CHL/2018. Accessed on 10/02/2019) or for the year 2016 for the west coast box. We use extract daily time series for each grid cell that falls within the NWA, SAB and CCS model regions, identify all pairs of consecutive days with valid data and convert changes in Chl to a carbon flux by using a fixed  C:Chl ratio of 60 (Ref63). The resulting weekly averages are interpolated to daily time series. The resulting change in DIC is significantly smaller than other sources and sinks in the current study. 

Winds are prescribed to 7 m/s in the summer and 10.5 m/s in the winter, to be consistent with representative wind data from the NOAA/Seawinds blended wind data set62.  Winter mixing in NWA is simulated by adding 0.1 mmol m-3 DIC daily from October through February (over 155 days). The value is based on a scaling analysis of vertical DIC gradients and diffusivity estimates in the region based on observational in 2015 in the NWA region and 2016 in the CCS. Literature values of physical vertical diffusion combined with vertical profiles of DIC from the ECOA 2015 cruse suggests that diffusive transports of DIC is negligible during summer. We use constant winds to minimize noise and to isolate the effect by changes in solubility on the carbonate system. Atmospheric pCO2 is set to 395 μatm. All carbonate equilibria and concentrations are calculated using the CO2SYS package[^29].

## Box model initiation and iteration
The model is initiated with the carbonate system defined by pCO2 in equilibrium with atmospheric values and alkalinity defined by initial salinity. The carbonate system is defined from DIC and TA and re-adjusted to changes in temperature and salinity at each time step using the CO2SYS package. pCO2_t(aq) is calculated using CO2SYS as well. ΔDICBio and ΔDICVertical are prescribed while ΔpCO2_Air-Sea is calculated at each time step using the relationship ΔpCO2_Air-Sea = 0.24 * kappa * K0 * (pCO2_t(aq) - (pCO2_(atm)). K0 is the CO2 gas solubility[^2] and kappa is defined as kappa = 0.251 * W2 * (Sc/660)-0.5 where W is wind speed in m/s and Sc is the Schmidt number[^65,^66].  Total alkalinity is calculated from prescribed salinity[^67].  The model is initiated to an equilibrium where pCO2_t0(aq) = pCO2_(atm), TA_t0 is based on salinity and temperature at t=0, and DIC is calculated from CO2SYS using pCO2(aq),  TAt0, temperature and salinity at t = 0. The model is spun up for 180 days using prescribed physical and biological conditions. The model iterated in time using the iteration

```
DICt+1 = DICt + ΔpCO2_Air-Sea + ΔDICBio + ΔDICVertical.
```

where each timestep is specifically conducted with the steps

```
1. Calculate the carbonate system defined by temperature, salinity, DIC, and TA at t_n.
2. Apply changes to the DIC concentration at t_n+1 by biological processes, 
vertical mixing and air-sea exchange (based on pCO2 values at t_n).
3. Recalculate the carbonate system using temperature, salinity, and TA at t_n,  and DIC at t_n+1.
4. Record pH, pCO2, and Ωarag at t_n+1.
```

The reason to use temperature, salinity, and TA at t_n  when performing step 3 is to make sure that all carbon parameters are calculated using the same physical conditions. We also calculate DICeq, a property that is based on the same TA, salinity, and temperature as the model but with pCO2(aq) relaxed to pCO2(atm) which can be interpreted as the air-sea exchange being infinitely fast. Here for simplicity, we simply assume pCO2(atm)  equals the dry CO2 mole fraction. 


## Installing / Getting started

Install the model with the following steps. The code and approaches has been tested on macOS and Linux and require a functioning python3.6 installed.

```shell
python3 -m venv /tmp/pco2box
source /tmp/pco2box/bin/activate
curl -O https://rsg.pml.ac.uk/shared_files/brj/pco2box.zip
unzip pco2box.zip
cd pco2box
pip install .
```

Run the model with

```python
import pco2box

pco2box.plot_perturbation(pert="ca") # California Current condtitions
pco2box.plot_perturbation(pert="ne") # North West Atlantic conditions
pco2box.plot_tmat_wind_mld()
```

This should create a folder called `pCO2box_figs` with the figures used in the paper.

## Links
- Publication where the model is used: https://doi.org/10.1038/s41467-020-16530-z

## Licensing
The code in this project is licensed under MIT license.
