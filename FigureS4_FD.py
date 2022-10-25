# Figure S4: Python script for base analysis and figures for Hartwick+2022
# Author: Victoria Hartwick
# Last Modified : 06/29/22

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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS4.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S4
# Wind turbines generate significant energy on an annual basis. Solar panel
# arrays would need to exceed 600 m2 on average to match energy produced by one
# turbine. (a) Annual energy production (GWh) by the Enercon E33 wind turbine
# over one Mars year compared with (b) the solar panel array dimensions (m2)
# needed to generate an equivalent amount of energy.
# -------------------------------------------------------------------------------
ct = 4
print('Making Figure S4')

# Plot
energylim=[0,15,24,35,100]  # energy requirement per sol (Rucker 2015 https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20150000526.pdf)
slevels = np.arange(0,1600,100)
smin = np.arange(0,1600,400)
fig, ax = py.subplots(1,2, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)

fig1 = ax[0].contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.aep_MY24[:,:]*1e-6, 20, cmap=py.cm.viridis)
ax[0].set_title('(a) AEP [GWh]')
ax[0].set_ylabel('Latitude')
ax[0].set_xlabel('Longitude')
ax[0].xaxis.set_major_locator(MultipleLocator(60))
ax[0].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig1,ax=ax[0])
cbar.set_label('[GWh]')

cmap = matplotlib.cm.plasma
norm = matplotlib.colors.BoundaryNorm(slevels,cmap.N,extend='max')
fig2 = ax[1].contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.S_solar[:,:], levels=slevels,cmap=py.cm.magma, norm=norm, extend='max')
ax[1].set_title('(b) Solar Array Size $\ [m^{2}$]')
cbar = py.colorbar(fig2,extend='max',ax=ax[1])
cbar.set_label('$\ [m^{2}$]')
fig2 = ax[1].contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.S_solar[:,:], levels=smin, cmap=py.cm.Greys)
py.clabel(fig2, fontsize=9,inline=1, fmt='%1.0f')
ax[1].set_xlabel('Longitude')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

