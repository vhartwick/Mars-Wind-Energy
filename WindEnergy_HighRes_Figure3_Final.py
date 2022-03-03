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

dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY28_highres/MY28_highres_savefile.nc'
DDS_MY28 = xr.open_dataset(dataDIR, decode_times=False)

# locate cardinal seasons in Ls array
areo = DDS_MY28.areo[:,0]
# note number of years
yrs = [0,1,2,3]


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
wpd_season_MY28 = []
swflx_season_MY28 = []
taudust_season_MY28 = []

#WPD
wpd_season_MY28 += [DDS_MY28.wpd[ls90,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls90')]
wpd_season_MY28 += [DDS_MY28.wpd[ls180,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls180')]
wpd_season_MY28 += [DDS_MY28.wpd[ls270,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls270')]
wpd_season_MY28 += [DDS_MY28.wpd[ls0,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls0')]
wpd_MY28 = xr.merge(wpd_season_MY28)

# SOLAR
swflx_season_MY28 += [DDS_MY28.swflx[ls90,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls90')]
swflx_season_MY28 += [DDS_MY28.swflx[ls180,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls180')]
swflx_season_MY28 += [DDS_MY28.swflx[ls270,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls270')]
swflx_season_MY28 += [DDS_MY28.swflx[ls0,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls0')]
swflx_MY28 = xr.merge(swflx_season_MY28)

# DUST OPTICAL DEPTH
taudust_season_MY28 += [DDS_MY28.taudust_VIS[ls90,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls90')]
taudust_season_MY28 += [DDS_MY28.taudust_VIS[ls180,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls180')]
taudust_season_MY28 += [DDS_MY28.taudust_VIS[ls270,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls270')]
taudust_season_MY28 += [DDS_MY28.taudust_VIS[ls0,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls0')]
taudust_MY28 = xr.merge(taudust_season_MY28)

# Same Thing for MY24
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
wpd_season_MY24 += [DDS_MY24.wpd[ls90,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls90')]
wpd_season_MY24 += [DDS_MY24.wpd[ls180,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls180')]
wpd_season_MY24 += [DDS_MY24.wpd[ls270,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls270')]
wpd_season_MY24 += [DDS_MY24.wpd[ls0,:,:,:,:].mean('time').mean('time_of_day_24').rename('wpd_Ls0')]
wpd_MY24 = xr.merge(wpd_season_MY24)

# SOLAR
swflx_season_MY24 += [DDS_MY24.swflx[ls90,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls90')]
swflx_season_MY24 += [DDS_MY24.swflx[ls180,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls180')]
swflx_season_MY24 += [DDS_MY24.swflx[ls270,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls270')]
swflx_season_MY24 += [DDS_MY24.swflx[ls0,:,:,:].mean('time').mean('time_of_day_24').rename('swflx_Ls0')]
swflx_MY24 = xr.merge(swflx_season_MY24)

# DUST OPTICAL DEPTH
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls90,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls90')]
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls180,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls180')]
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls270,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls270')]
taudust_season_MY24 += [DDS_MY24.taudust_VIS[ls0,:,:,:].mean('time').mean('time_of_day_24').rename('taudust_Ls0')]
taudust_MY24 = xr.merge(taudust_season_MY24)

# -------------------------------------------------------------------------------
# FIGURE 3
# Wind energy increases during global dust storms while solar power is reduced. 
# Difference in (A) wind power density, (B) solar power density, and (C) dust 
# optical depth for a global storm year (MY 28) minus a non-global storm year 
# (MY24) at perihelion (Ls=260-280). 
# -------------------------------------------------------------------------------
ct = 3
print('Making Figure 3')

# just plot dust storm Ls
levels_wpd_diff = np.arange(-120,130,10)
levels_spd_diff = [-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60]
levels_dod_diff = [-1.6,-1.4,-1.2,-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1.,1.2,1.4,1.6]
levels_dod_diff = np.arange(-5,5.25,0.25)

# panel a : wind profile density difference
fig, ([ax1,ax2,ax3]) = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)
fig1 = ax1.contourf(DDS_MY24.lon, DDS_MY24.lat, wpd_MY28.wpd_Ls270.sel(zagl=50)-
                    wpd_MY24.wpd_Ls270.sel(zagl=50),cmap=py.cm.RdBu_r)
ax1.set_title('(A) WPD')
ax1.set_ylabel('Latitude')
ax1.set_xlabel('Longitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig1, ax=ax1,orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')
print(max(wpd_MY28.wpd_Ls270.sel(zagl=50).values-wpd_MY24.wpd_Ls270.sel(zagl=50).values))
print(min(swflx_MY28.swflx_Ls270.values-swflx_MY24.swflx_Ls270.values))

# panel b - solar power density difference
fig2 = ax2.contourf(DDS_MY24.lon, DDS_MY24.lat, swflx_MY28.swflx_Ls270-swflx_MY24.swflx_Ls270, 
                    levels=levels_spd_diff,cmap=py.cm.RdBu_r)
ax2.set_title('(B) SPD')
ax2.set_xlabel('Longitude')
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig2, ax=ax2,orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel c - dust optical depth difference
fig3 = ax3.contourf(DDS_MY24.lon, DDS_MY24.lat, taudust_MY28.taudust_Ls270-taudust_MY24.taudust_Ls270, 
                    levels=levels_dod_diff,cmap=py.cm.RdBu_r)
ax3.set_title('(C) Dust Optical Depth')
ax3.set_xlabel('Longitude')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
py.colorbar(fig3, ax=ax3,orientation='horizontal', ticks=[-5,-2.5,0,2.5,5])

py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300) 

