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
dataDIR = '/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_solarwindpercent300.nc'
DDS_300 = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_windpercent2.nc'
DDS_PCT = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_savefile.nc'
DDS_swflx = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S7
# (As in Figure 6a-c) Percent time of MY24 in hourly intervals that wind power, 
# solar power and wind and solar combined can generates >=24kW. Solar power is 
# calculated based on a 300 m2 array (29) with Cp=0.2.
# -------------------------------------------------------------------------------
ct = 7
print('Making Figure S7')

time = len(DDS_swflx.time)*len(DDS_swflx.time_of_day_24)
levels = np.arange(0,100,5)
level_lim = [0,50,1000]


fig, axs = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5)) # define number of rows, columns
fig.subplots_adjust(wspace=0.05, hspace=0.01)

# panel a - wind
cs = axs[0].contourf(DDS_PCT.lon, DDS_PCT.lat, (DDS_PCT.count_wind24/time)*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[0].contour(DDS_PCT.lon, DDS_PCT.lat, (DDS_PCT.count_wind24/time)*100, colors='k', levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[0].yaxis.set_major_locator(MultipleLocator(30))
axs[0].set_title('(A) Wind')
axs[0].set_xlabel('Longitude')
axs[0].set_ylabel('Latitude')


# panel b - solar
cs = axs[1].contourf(DDS_300.lon, DDS_300.lat, DDS_300.percent_solar300*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[1].contour(DDS_300.lon, DDS_300.lat, DDS_300.percent_solar300*100, levels=levels, colors='k')
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[1].set_title('(B) Solar, $\ 300m^{2}$')
axs[1].set_xlabel('Longitude')

# panel c - total
cs = axs[2].contourf(DDS_300.lon, DDS_300.lat, DDS_300.percent_total300*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[2].contour(DDS_300.lon, DDS_300.lat, DDS_300.percent_total300*100, colors='k',levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[2].set_title('(C) Solar + Wind')
axs[2].set_xlabel('Longitude')

norm= matplotlib.colors.Normalize(vmin=cs.cvalues.min(), vmax=cs.cvalues.max())
sm = py.cm.ScalarMappable(norm=norm, cmap = cs.cmap)
sm.set_array([])
cbar_ax = py.gcf().add_axes([0.92, 0.125, 0.01, 0.755]) #left, bottom, width, height
clb = fig.colorbar(sm, cax=cbar_ax, orientation='vertical', label='% year')


py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300) 
