# Figure1 : Python script for base analysis and figures for Hartwick+2022
# Author: Victoria Hartwick
# Last Modified : 6/29/22

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
dataDIR = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/Figure1.nc'
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

levels_swf = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
levels_zsurf=np.linspace(min(DDS_MY24.zsurf),max(DDS_MY24.zsurf),5)
print(levels_zsurf)
fig, axs = py.subplots(2,6, sharey=True, figsize=(12,10))
fig.subplots_adjust(wspace=0.1,hspace=0.2)

gs = gridspec.GridSpec(2,6)
ax1 = py.subplot(gs[0,:3])
ax2 = py.subplot(gs[0,3:])
ax3 = py.subplot(gs[1,0:2])
ax4 = py.subplot(gs[1,2:4])
ax5 = py.subplot(gs[1,4:])

print(DDS_MY24.zsurf)
# panel a - 50m, consistent with base turbine hub height
fig1 = ax1.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_annual_ave_MY24.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax1.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,levels=levels_zsurf, cmap=py.cm.Greys_r)
ax1.set_title('(a) 50m')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig1,ax=ax1)
clb.set_label('$\ [W/m^{2}$]')

# panel b - wind (50m) to solar ratio
fig2 = ax2.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_annual_ave_MY24.sel(zagl=50)/DDS_MY24.swflx_annual_ave_MY24, cmap=py.cm.viridis)
ax2.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,levels=levels_zsurf, cmap=py.cm.Greys_r)
ax2.set_title('(b) 50m Wind:Solar')
ax2.set_xlabel('Longitude')
ax2.set_yticklabels([])
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig2, ax=ax2)

# panel c - 5m
fig3 = ax3.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_annual_ave_MY24.sel(zagl=5), cmap=py.cm.viridis, levels=levels_swf)
ax3.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,levels=levels_zsurf, cmap=py.cm.Greys_r)
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')
ax3.set_title('(c) 5m')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig3, ax=ax3, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel d - 30m
fig4 = ax4.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_annual_ave_MY24.sel(zagl=30), cmap=py.cm.viridis, levels=levels_swf)
ax4.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,levels=levels_zsurf, cmap=py.cm.Greys_r)
ax4.set_xlabel('Longitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
ax4.set_title('(d) 30m')
ax4.set_yticklabels([])
clb=py.colorbar(fig4, ax=ax4, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel e - 100m
fig5 = ax5.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_annual_ave_MY24.sel(zagl=100), cmap=py.cm.viridis, levels=levels_swf)
ax5.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf,levels=levels_zsurf, cmap=py.cm.Greys_r)
ax5.set_xlabel('Longitude')
ax5.set_title('(e) 100m')
ax5.xaxis.set_major_locator(MultipleLocator(60))
ax5.yaxis.set_major_locator(MultipleLocator(30))
ax5.set_yticklabels([])
clb = py.colorbar(fig5, ax=ax5, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300)

