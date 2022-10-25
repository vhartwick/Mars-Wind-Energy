# Figur4 : Python script for base analysis and figures for Hartwick+2022
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
dataDIR = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Figure4.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE 4
# Wind power is greatest at night when solar power is at a minimum. Panels show
# the annual average (A) daytime (7.5-16.5LT) and (B) nighttime wind power
# density and (C) their ratio.
# -------------------------------------------------------------------------------
ct = 4
print('Making Figure 4')

levels=np.arange(0,220,20)
level_ratio = np.arange(0,5.5,0.5)
level_ratio = np.arange(0,10.5,0.5)

fig, axs = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)
axs = axs.ravel()

# panel a - nighttime winds
fig1 = axs[0].contourf(DDS_MY24.lon,DDS_MY24.lat,DDS_MY24.wpd_night,levels=levels,cmap=py.cm.viridis)
axs[0].set_title('(a) Nighttime Wind')
axs[0].set_xlabel('Longitude')
axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[0].yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig1,ax=axs[0],orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')
axs[0].set_ylabel('Latitude')

# panel b - daytime winds
fig = axs[1].contourf(DDS_MY24.lon,DDS_MY24.lat,DDS_MY24.wpd_day,levels=levels,cmap=py.cm.viridis)
axs[1].set_title('(b) Daytime Wind')
axs[1].set_xlabel('Longitude')
clb=py.colorbar(fig,ax=axs[1],orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel c - difference
cmap = matplotlib.cm.plasma
norm = matplotlib.colors.BoundaryNorm(level_ratio,cmap.N,extend='max')

fig = axs[2].contourf(DDS_MY24.lon,DDS_MY24.lat,DDS_MY24.wpd_night/DDS_MY24.wpd_day,cmap=py.cm.plasma,levels=level_ratio,extend='max',norm=norm)
axs[2].set_title('(c) Night:Day')
axs[2].set_xlabel('Longitude')
py.colorbar(fig,ax=axs[2],orientation='horizontal',extend='max')

py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300)
