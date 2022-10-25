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
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_savefile.nc'
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

# note number of years
yrs = [0,1,2,3]

# locate cardinal seasons in Ls array
areo = DDS_MY24.areo[:,0]
lsrange, lsrange1, lsrange2, lsrange3 = [], [], [], []
for i in yrs:
    # classic Ls ranges
    loc = np.where(np.logical_and(areo >= 360*i + 80, areo <= 360*i + 100.))
    loc1 = np.where(np.logical_and(areo >= 360*i + 170, areo <= 360*i + 190.))
    loc2 = np.where(np.logical_and(areo >= 360*i + 260, areo <= 360*i + 280.))
    loc3 = np.where(np.logical_and(areo >= 360*i + 350, areo <= 360*(i+1) + 10.))

    lsrange.append(loc[0])
    lsrange1.append(loc1[0])
    lsrange2.append(loc2[0])
    lsrange3.append(loc3[0])

ls90 = np.concatenate(lsrange)
ls180 = np.concatenate(lsrange1)
ls270 = np.concatenate(lsrange2)
ls0 = np.concatenate(lsrange3)

# average all years to produce multiyear seasonal averages (MC_season)
wpd_season_MY24 = []
swflx_season_MY24 = []
taudust_season_MY24 = []

#WPD
wpd_season_MY24 += [DDS_MY24.wpd[ls90,:,:,:,:].mean(('time','time_of_day_24')).rename('wpd_Ls90')]
wpd_season_MY24 += [DDS_MY24.wpd[ls180,:,:,:,:].mean(('time','time_of_day_24')).rename('wpd_Ls180')]
wpd_season_MY24 += [DDS_MY24.wpd[ls270,:,:,:,:].mean(('time','time_of_day_24')).rename('wpd_Ls270')]
wpd_season_MY24 += [DDS_MY24.wpd[ls0,:,:,:,:].mean(('time','time_of_day_24')).rename('wpd_Ls0')]
wpd_MY24 = xr.merge(wpd_season_MY24)

# SOLAR
swflx_season_MY24 += [DDS_MY24.swflx[ls90,:,:,:].mean(('time','time_of_day_24')).rename('swflx_Ls90')]
swflx_season_MY24 += [DDS_MY24.swflx[ls180,:,:,:].mean(('time','time_of_day_24')).rename('swflx_Ls180')]
swflx_season_MY24 += [DDS_MY24.swflx[ls270,:,:,:].mean(('time','time_of_day_24')).rename('swflx_Ls270')]
swflx_season_MY24 += [DDS_MY24.swflx[ls0,:,:,:].mean(('time','time_of_day_24')).rename('swflx_Ls0')]
swflx_MY24 = xr.merge(swflx_season_MY24)

# DUST OPTICAL DEPTH
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls90,:,:,:].mean(('time','time_of_day_24')).rename('taudust_Ls90')]
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls180,:,:,:].mean(('time','time_of_day_24')).rename('taudust_Ls180')]
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls270,:,:,:].mean(('time','time_of_day_24')).rename('taudust_Ls270')]
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls0,:,:,:].mean(('time','time_of_day_24')).rename('taudust_Ls0')]
taudust_MY24 = xr.merge(taudust_season_MY24)

# Cardinal Seasons FV3 (Using DIURN FILE)
import warnings; warnings.simplefilter('ignore')

levels_swf = np.arange(0,160,10)
levels_swf2 = np.arange(0,240,20)
level_lim = [100,500]
dust_lim = [0,0.5] #is there a way to calculate this? optical depth when solar cell starts working at half efficiency or something?

# Plot WPD
fig, ([ax1,ax2],[ax3,ax4]) = py.subplots(2,2, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)
#fig.suptitle('Seasonal Average Wind Power Density at 50m $\ [W/m^{2}]$, MY 24')

# panel a - Ls90
fig1 = ax1.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_MY24.wpd_Ls90.sel(zagl=50).squeeze(), levels=levels_swf,cmap=py.cm.viridis)
ax1.set_title('(a) $\ L_{s}=80-100\degree$, NH Summer Solstice')
ax1.set_ylabel('Latitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig1, ax=ax1)
clb.set_label('$\ [W/m^{2}$]')

fig1 = ax1.contour(DDS_MY24.lon, DDS_MY24.lat, swflx_MY24.swflx_Ls90, levels=levels_swf2, cmap=py.cm.plasma)
clabel = ax1.clabel(fig1, fontsize=9,inline=1, fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]

# panel b - Ls180
fig2 = ax2.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_MY24.wpd_Ls180.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax2.set_title('$\ (b) L_{s}=170-190\degree$, NH Fall Equinox')
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
clb= py.colorbar(fig2, ax=ax2)
clb.set_label('$\ [W/m^{2}$]')

fig2 = ax2.contour(DDS_MY24.lon, DDS_MY24.lat, swflx_MY24.swflx_Ls180, levels=levels_swf2, cmap=py.cm.plasma)
clabel = ax2.clabel(fig2, fontsize=9,fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]

# panel c - Ls270
fig3 = ax3.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_MY24.wpd_Ls270.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax3.set_title('$\ (c) L_{s}=260-280\degree$, NH Winter Solstice')
ax3.set_ylabel('Latitude')
ax3.set_xlabel('Longitude')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig3, ax=ax3)
clb.set_label('$\ [W/m^{2}$]')

fig3 = ax3.contour(DDS_MY24.lon, DDS_MY24.lat, swflx_MY24.swflx_Ls270, levels=levels_swf2, cmap=py.cm.plasma)
clabel=ax3.clabel(fig3, fontsize=9,inline=1, fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]

# panel d - Ls0
fig4 = ax4.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_MY24.wpd_Ls0.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax4.set_title('$\ (d) L_{s}=350-10\degree$, NH Spring Equinox')
ax4.set_xlabel('Longitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig4, ax=ax4)
clb.set_label('$\ [W/m^{2}$]')

fig4 = ax4.contour(DDS_MY24.lon, DDS_MY24.lat, swflx_MY24.swflx_Ls0, levels=levels_swf2, cmap=py.cm.plasma)
clabel=ax4.clabel(fig4, fontsize=9,inline=1, fmt='%1.0f')
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabel]


py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300)

# Send out Data
var1 = xr.DataArray(wpd_MY24.wpd_Ls90.sel(zagl=50), coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var1.name = 'wpd_Ls90'

var2 = xr.DataArray(wpd_MY24.wpd_Ls180.sel(zagl=50), coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var2.name = 'wpd_Ls180'

var3 = xr.DataArray(wpd_MY24.wpd_Ls270.sel(zagl=50), coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var3.name = 'wpd_Ls270'

var4 = xr.DataArray(wpd_MY24.wpd_Ls0.sel(zagl=50), coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var4.name = 'wpd_Ls0'

var5 = xr.DataArray(swflx_MY24.swflx_Ls90, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var5.name = 'swflx_Ls90'

var6 = xr.DataArray(swflx_MY24.swflx_Ls180, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var6.name = 'swflx_Ls180'

var7 = xr.DataArray(swflx_MY24.swflx_Ls270, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var7.name = 'swflx_Ls270'

var8 = xr.DataArray(swflx_MY24.swflx_Ls0, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var8.name = 'swflx_Ls0'

NEW_DF = xr.merge([var1,var2,var3,var4,var5,var6,var7,var8])
NEW_DF.to_netcdf(f'{PATH}/Data/Figure{ct}.nc')


