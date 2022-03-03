# WindEnergy_Figures : Python script for base analysis and figures for Hartwick+2021a
# Author: Victoria Hartwick
# Created : 12/31/21

# Import Scientific Modules
import os
import matplotlib
import matplotlib.pyplot as py
import numpy as np
from matplotlib.ticker import MultipleLocator, FuncFormatter #format ticks
from numpy import sqrt, exp, max, mean, min, log, log10,int
import xarray as xr
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

# Get Path
PATH = os.getcwd()

# Import Base Datasets
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_savefile.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE 4
# Wind power is greatest at night when solar power is at a minimum. Panels show
# the annual average (A) daytime (7.5-16.5LT) and (B) nighttime wind power
# density and (C) their ratio. 
# -------------------------------------------------------------------------------
ct = 4
print('Making Figure 4')

# generate day and night masks
tmp = DDS_MY24.wpd.sel(zagl=50,time_of_day_24=slice(0,8))
tmp2 = DDS_MY24.wpd.sel(zagl=50,time_of_day_24=slice(18,24))
night = xr.concat([tmp,tmp2],dim='time_of_day_24')
night = night.mean(('time','time_of_day_24'))

day = DDS_MY24.wpd.sel(zagl=50,time_of_day_24=slice(8,17)).mean(('time','time_of_day_24'))
weights = np.cos(np.deg2rad(DDS_MY24.lat))
weights.name="weights"
ratio = night/day
ga = ratio.weighted(weights).mean(('lon','lat'))
print(ga)
# Plot
levels=np.arange(0,220,20)
level_ratio = np.arange(0,5.5,0.5)
level_ratio = np.arange(0,10.5,0.5)

fig, axs = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)
axs = axs.ravel()

# panel a - nighttime winds
fig1 = axs[0].contourf(DDS_MY24.lon,DDS_MY24.lat,night,levels=levels,cmap=py.cm.viridis)
axs[0].set_title('(A) Nighttime Wind')
axs[0].set_xlabel('Longitude')
axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[0].yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig1,ax=axs[0],orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')
axs[0].set_ylabel('Latitude')

# panel b - daytime winds
fig = axs[1].contourf(DDS_MY24.lon,DDS_MY24.lat,day,levels=levels,cmap=py.cm.viridis)
axs[1].set_title('(B) Daytime Wind')
axs[1].set_xlabel('Longitude')
clb=py.colorbar(fig,ax=axs[1],orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel c - difference
cmap = matplotlib.cm.plasma
norm = matplotlib.colors.BoundaryNorm(level_ratio,cmap.N,extend='max')

fig = axs[2].contourf(DDS_MY24.lon,DDS_MY24.lat,night/day,cmap=py.cm.plasma,levels=level_ratio,extend='max',norm=norm)
axs[2].set_title('(C) Night/Day')
axs[2].set_xlabel('Longitude')
py.colorbar(fig,ax=axs[2],orientation='horizontal',extend='max')

py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300) 

