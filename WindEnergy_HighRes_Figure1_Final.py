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
# FIGURE 1
# Wind power potential varies with altitude above the surface. At 50m the ratio 
# of wind to solar power can exceed one. Annual average (MY24) wind power density
# [W/m2] at (A) 5m, (B) 30m, (C) 100 m and (D) the reference altitude for our 
# base turbine, 50m above the surface. Panel (E) shows the ratio of wind to solar
# power at 50m. Solid contours in all panels show surface topography.
# -------------------------------------------------------------------------------
ct = 1
print('Making Figure 1')

# Find Annual Average
# First simulation year (MY24)
wpd_annual_ave_MY24 = DDS_MY24.wpd.mean('time').mean('time_of_day_24')
swflx_annual_ave_MY24 = DDS_MY24.swflx.mean('time').mean('time_of_day_24')
taudust_annual_ave_MY24 = DDS_MY24.taudust_VIS.mean('time').mean('time_of_day_24')

wpd_time_ave_MY24 = DDS_MY24.wpd.mean('time_of_day_24')
wpd_diurn_ave_MY24 = DDS_MY24.wpd.mean('time')
swflx_diurn_ave_MY24 = DDS_MY24.swflx.mean('time')
print(max(wpd_time_ave_MY24).values, max(wpd_diurn_ave_MY24).values, max(DDS_MY24.wpd).values)

levels_swf = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
fig, axs = py.subplots(2,6, sharey=True, figsize=(12,10))
fig.subplots_adjust(wspace=0.1,hspace=0.1)

gs = gridspec.GridSpec(2,6)
ax1 = py.subplot(gs[0,0:2])
ax2 = py.subplot(gs[0,2:4])
ax3 = py.subplot(gs[0,4:])
ax4 = py.subplot(gs[1,:3])
ax5 = py.subplot(gs[1,3:])

# panel a - 5m
fig1 = ax1.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_annual_ave_MY24.sel(zagl=5), cmap=py.cm.viridis, levels=levels_swf)
ax1.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,5, cmap=py.cm.Greys_r)
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
ax1.set_title('(A) 5m')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig1, ax=ax1, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel b - 30m
fig2 = ax2.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_annual_ave_MY24.sel(zagl=30), cmap=py.cm.viridis, levels=levels_swf)
ax2.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,5, cmap=py.cm.Greys_r)
ax2.set_xlabel('Longitude')
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
ax2.set_title('(B) 30m')
ax2.set_yticklabels([])
clb=py.colorbar(fig2, ax=ax2, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel c - 100m
fig3 = ax3.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_annual_ave_MY24.sel(zagl=100), cmap=py.cm.viridis, levels=levels_swf)
ax3.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,5, cmap=py.cm.Greys_r)
ax3.set_xlabel('Longitude')
ax3.set_title('(C) 100m')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
ax3.set_yticklabels([])
clb = py.colorbar(fig3, ax=ax3, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel d - 50m, consistent with base turbine hub height
fig4 = ax4.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_annual_ave_MY24.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax4.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,5, cmap=py.cm.Greys_r)
ax4.set_title('(D) 50m')
ax4.set_xlabel('Longitude')
ax4.set_ylabel('Latitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig4,ax=ax4)
clb.set_label('$\ [W/m^{2}$]')

# panel e - wind (50m) to solar ratio
fig5 = ax5.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_annual_ave_MY24.sel(zagl=50)/swflx_annual_ave_MY24, cmap=py.cm.viridis)
ax5.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,5, cmap=py.cm.Greys_r)
ax5.set_title('(E) 50m Wind:Solar')
ax5.set_xlabel('Longitude')
ax5.set_yticklabels([])
ax5.xaxis.set_major_locator(MultipleLocator(60))
ax5.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig5, ax=ax5)

py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300) 

