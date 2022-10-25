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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS1.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S1
# Zonal average solar power density (blue dashed line) and average wind power
# density in 1° longitude bands (colored lines) vs Ls at (a) 70N, (b) 50N,
# (c) 35N, (d) 43S and (e) 50S. Surface short wave flux [W/m2] is calculated in
# the model based on the atmospheric composition and scattering and absorption
# by atmospheric aerosols. Solar power is greatest in the summer hemisphere and
# is reduced during large scale dust-events (as at Ls=240 in panels d,e). Wind
# power regularly exceeds solar power in the winter hemisphere and therefore is
# well-suited to complement solar arrays to power future human missions.
# -------------------------------------------------------------------------------

ct = 1
print('Making Figure S1')

areo_MY24 = DDS_MY24.areo[:,0]-360*3

fig, ax = py.subplots(5,1, sharex=True, figsize=(12,10))
fig.subplots_adjust(wspace=0.05)

# panel a
fig1 = ax[0].plot(DDS_MY24.areo_MY24, DDS_MY24.diurn_wpd_MY24.sel(lat=70, zagl=50, method='nearest'))
fig1 = ax[0].plot(DDS_MY24.areo_MY24, DDS_MY24.zonal_diurn_swflx_MY24.sel(lat=70,method='nearest'), linewidth=3,linestyle='--')
ax[0].set_title('$\ (a) 70\degree$N')
ax[0].set_ylabel('$\ [W/m^{2}]$')
ax[0].xaxis.set_major_locator(MultipleLocator(60))
ax[0].set_xlim([0,360])

# panel b
fig1 = ax[1].plot(DDS_MY24.areo_MY24, DDS_MY24.diurn_wpd_MY24.sel(lat=50, zagl=50, method='nearest'))
fig1 = ax[1].plot(DDS_MY24.areo_MY24, DDS_MY24.zonal_diurn_swflx_MY24.sel(lat=50,method='nearest'),linewidth=3,linestyle='--')
ax[1].set_title('$\ (b) 50\degree$N')
ax[1].set_ylabel('$\ [W/m^{2}]$')

# panel c
fig1 = ax[2].plot(DDS_MY24.areo_MY24, DDS_MY24.diurn_wpd_MY24.sel(lat=35, zagl=50, method='nearest'))
fig1 = ax[2].plot(DDS_MY24.areo_MY24, DDS_MY24.zonal_diurn_swflx_MY24.sel(lat=35,method='nearest'),linewidth=3,linestyle='--')
ax[2].set_title('$\ (c) 35\degree$N')
ax[2].set_ylabel('$\ [W/m^{2}]$')

# panel d
fig2 = ax[3].plot(DDS_MY24.areo_MY24, DDS_MY24.diurn_wpd_MY24.sel(lat=-43, zagl=50, method='nearest'))
fig2 = ax[3].plot(DDS_MY24.areo_MY24, DDS_MY24.zonal_diurn_swflx_MY24.sel(lat=-43,  method='nearest'),linewidth=3,linestyle='--')
ax[3].set_title('$\ (d) 43\degree$S')
ax[3].set_ylabel('$\ [W/m^{2}]$')

# panel e
#a = np.mean(tmp_SY1[:,:,:,24:33],axis=3)   # average between lon = 90-120 E
fig2 = ax[4].plot(DDS_MY24.areo_MY24, DDS_MY24.diurn_wpd_MY24.sel(lat=-50, zagl=50, method='nearest'))
fig2 = ax[4].plot(DDS_MY24.areo_MY24, DDS_MY24.zonal_diurn_swflx_MY24.sel(lat=-50,  method='nearest'),linewidth=3,linestyle='--')
ax[4].set_title('$\ (e) 50\degree$S')
ax[4].set_ylabel('$\ [W/m^{2}]$')
ax[4].set_xlabel('$\ L_{s}$')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

