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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS12.nc'
DDS= xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S12
# -------------------------------------------------------------------------------
ct = 12
print('Making Figure S12')

fig, ax = py.subplots(1,2,figsize=(12,5))
fig.subplots_adjust(wspace=0.15, hspace = 0.15)
ax = ax.ravel()
.84,0,0,0,0,0,0]

# panel a - Opportunity
fig1 = ax[0].plot(DDS.tod,DDS.opppow, color='purple')
ax[0].set_xlabel('Time of Day (LT)')
ax[0].set_ylabel('[W]')
ax[0].axhline(1, color='red', linestyle='--', label='Energy Requirement')
ax[0].set_title('(a) Aeolos-V Power Output [W], $\ L_{s}=190\degree$')

# panel b - AEP
fig2 = ax[1].contourf(DDS.lon, DDS.lat, DDS.aep, cmap=py.cm.viridis)
ax[1].set_ylabel('Latitude')
ax[1].set_xlabel('Longitude')
ax[1].set_title('(b) Aeolos-V AEP [kWh]')
ax[1].xaxis.set_major_locator(MultipleLocator(60))
ax[1].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig2,ax=ax[1])
cbar.set_label('[kWh]')
py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}_alt.eps',dpi=300)

