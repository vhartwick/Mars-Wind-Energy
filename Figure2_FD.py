#ndEnergy_Figures : Python script for base analysis and figures for Hartwick+2022
# Author: Victoria Hartwick
# Last Modified : 6/29/22

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
dataDIR = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Figure2.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE 2
# Wind power density exceeds solar power density particularly in the winter
# hemisphere mid- to polar latitudes. Seasonal wind power density [W/m2] at 50m
# for the Mars cardinal seasons in 20° Ls bins. Solid contours show the surface
# solar power density for the same season. Solar power varies latitudinally with
# season. Baroclinic wave activity adjacent to the winter polar vortex accelerates
# atmospheric wind speeds when solar energy yields are depleted.
# -------------------------------------------------------------------------------
ct = 2
print('Making Figure 2')
levels_swf = np.arange(0,160,10)
levels_swf2 = np.arange(0,240,20)
level_lim = [100,500]
dust_lim = [0,0.5] #is there a way to calculate this? optical depth when solar cell starts working at half efficiency or something?

# Plot WPD
fig, ([ax1,ax2],[ax3,ax4]) = py.subplots(2,2, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)
#fig.suptitle('Seasonal Average Wind Power Density at 50m $\ [W/m^{2}]$, MY 24')

# panel a - Ls90
fig1 = ax1.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_Ls90.sel(zagl=50).squeeze(), levels=levels_swf,cmap=py.cm.viridis)
ax1.set_title('(a) $\ L_{s}=80-100\degree$, NH Summer Solstice')
ax1.set_ylabel('Latitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig1, ax=ax1)
clb.set_label('$\ [W/m^{2}$]')

fig1 = ax1.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.swflx_Ls90, levels=levels_swf2, cmap=py.cm.plasma)
clabel = ax1.clabel(fig1, fontsize=9,inline=1, fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]

# panel b - Ls180
fig2 = ax2.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_Ls180.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax2.set_title('$\ (b) L_{s}=170-190\degree$, NH Fall Equinox')
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
clb= py.colorbar(fig2, ax=ax2)
clb.set_label('$\ [W/m^{2}$]')

fig2 = ax2.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.swflx_Ls180, levels=levels_swf2, cmap=py.cm.plasma)
clabel = ax2.clabel(fig2, fontsize=9,fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]

# panel c - Ls270
fig3 = ax3.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_Ls270.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax3.set_title('$\ (c) L_{s}=260-280\degree$, NH Winter Solstice')
ax3.set_ylabel('Latitude')
ax3.set_xlabel('Longitude')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig3, ax=ax3)
clb.set_label('$\ [W/m^{2}$]')

fig3 = ax3.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.swflx_Ls270, levels=levels_swf2, cmap=py.cm.plasma)
clabel=ax3.clabel(fig3, fontsize=9,inline=1, fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]

# panel d - Ls0
fig4 = ax4.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.wpd_Ls0.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax4.set_title('$\ (d) L_{s}=350-10\degree$, NH Spring Equinox')
ax4.set_xlabel('Longitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig4, ax=ax4)
clb.set_label('$\ [W/m^{2}$]')

fig4 = ax4.contour(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.swflx_Ls0, levels=levels_swf2, cmap=py.cm.plasma)
clabel=ax4.clabel(fig4, fontsize=9,inline=1, fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]


py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300)

