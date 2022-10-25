# Figure S7: Python script for base analysis and figures for Hartwick+2022
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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS7.nc'
DDS_PCT = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S7
# (As in Figure 6a-c) Percent time of MY24 in hourly intervals that wind power,
# solar power and wind and solar combined can generates >=24kW. Solar power is
# calculated based on a 300 m2 array (29) with Cp=0.2.
# -------------------------------------------------------------------------------
ct = 7
print('Making Figure S7')

levels = np.arange(0,100,5)
level_lim = [0,50,1000]

fig, axs = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5)) # define number of rows, columns
fig.subplots_adjust(wspace=0.05, hspace=0.01)

# panel a - wind
cs = axs[0].contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_wind24*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[0].contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_wind24*100, colors='k', levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[0].yaxis.set_major_locator(MultipleLocator(30))
axs[0].set_title('(a) Wind')
axs[0].set_xlabel('Longitude')
axs[0].set_ylabel('Latitude')


# panel b - solar
cs = axs[1].contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_solar300*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[1].contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_solar300*100, levels=levels, colors='k')
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[1].set_title('(b) Solar, $\ 300m^{2}$')
axs[1].set_xlabel('Longitude')

# panel c - total
cs = axs[2].contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_total300*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[2].contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_total300*100, colors='k',levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[2].set_title('(c) Solar + Wind')
axs[2].set_xlabel('Longitude')

norm= matplotlib.colors.Normalize(vmin=cs.cvalues.min(), vmax=cs.cvalues.max())
sm = py.cm.ScalarMappable(norm=norm, cmap = cs.cmap)
sm.set_array([])
cbar_ax = py.gcf().add_axes([0.92, 0.125, 0.01, 0.755]) #left, bottom, width, height
clb = fig.colorbar(sm, cax=cbar_ax, orientation='vertical', label='% year')


py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

