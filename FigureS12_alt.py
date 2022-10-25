# Figure S12: Python script for base analysis and figures for Hartwick+2022
# Author: Victoria Hartwick
# Last Modified: 6/29/22

# Import Scientific Modules
import os
import matplotlib
import matplotlib.pyplot as py
import numpy as np
from matplotlib.ticker import MultipleLocator, FuncFormatter #format ticks
from numpy import sqrt, exp, max, mean, min, log, log10
import xarray as xr
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

# Get Path
PATH = os.getcwd()

# Import Base Datasets
dataDIR = '/lou/la3/mkahre/MCMC/tmp/OLIVIA/AAV_300W_power_MY24_HIGHRES.nc' 
DDS= xr.open_dataset(dataDIR, decode_times=False)

# locate indeces of ROI lon/lat

lon_roi = [2.28]
lat_roi = [-5.23]

lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS.lon-tmp[i])).argmin()
    lonlim_idx += [idx] 

tmp = np.array(lat_roi) 
for i in range(0,len(tmp)):
    idx = (np.abs(DDS.lat-tmp[i])).argmin()
    latlim_idx += [idx] 
    
# convert kW to W by multiplying by 1000
ls190_exact = DDS.high_res_pow_diurn_SY1[77,:,22,0] * 1000
print(ls190_exact)
# -------------------------------------------------------------------------------
# FIGURE S12
# -------------------------------------------------------------------------------
ct = 12
print('Making Figure S12')

fig, ax = py.subplots(1,2,figsize=(12,5))
fig.subplots_adjust(wspace=0.15, hspace = 0.15)
ax = ax.ravel()

tod = np.arange(0.5,24.5,1)
print(tod)
opppow = [0, 3.02, 4.79, 6.92, 13.92, 13.92, 25.1,43.02,43.02,13.92,2.84,0,0,4.88,6.92,13.84,6.92,2.84,0,0,0,0,0,0]

# panel a - Opportunity
fig1 = ax[0].plot(tod,opppow, color='purple')
ax[0].set_xlabel('Time of Day (LT)')
ax[0].set_ylabel('[W]')
ax[0].axhline(1, color='red', linestyle='--', label='Energy Requirement')
ax[0].set_title('(a) Aeolos-V Power Output [W], $\ L_{s}=190\degree$')

# panel b - AEP
aep = np.sum(DDS.high_res_pow_diurn_SY1, axis=0)
aep = np.sum(aep, axis=0)*(5/1000)

# Calculate global average AEP
weights = np.cos(np.deg2rad(DDS.lat))
weights.name = "weights"
print('global average aep=', aep.weighted(weights).mean(('lon','lat')))



fig2 = ax[1].contourf(DDS.lon, DDS.lat, aep, cmap=py.cm.viridis)
ax[1].set_ylabel('Latitude')
ax[1].set_xlabel('Longitude')
ax[1].set_title('(b) Aeolos-V AEP [kWh]')
ax[1].xaxis.set_major_locator(MultipleLocator(60))
ax[1].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig2,ax=ax[1])
cbar.set_label('[kWh]')
py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}_alt.eps',dpi=300)

# Save Data
var1 = xr.DataArray(aep, coords={'lat':DDS.lat, 'lon':DDS.lon},dims=['lat','lon'])

NEW_DF = xr.Dataset({'tod':tod,'opppow':opppow, 'aep':var1})
NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')
