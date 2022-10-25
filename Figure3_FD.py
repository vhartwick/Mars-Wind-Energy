# Figure3 : Python script for base analysis and figures for Hartwick+2022
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
dataDIR = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Figure3.py'
DDS = xr.open_dataset(dataDIR, decode_times=False)

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
fig1 = ax1.contourf(DDS.lon, DDS.lat, DDS.wpd_MY28_Ls270.sel(zagl=50)-DDS.wpd_MY24_Ls270.sel(zagl=50),
                    cmap=py.cm.RdBu_r)
ax1.set_title('(a) WPD')
ax1.set_ylabel('Latitude')
ax1.set_xlabel('Longitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig1, ax=ax1,orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel b - solar power density difference
fig2 = ax2.contourf(DDS.lon, DDS.lat, DDS.swflx_MY28_Ls270-DDS.swflx_MY24_Ls270,
                    levels=levels_spd_diff,cmap=py.cm.RdBu_r)
ax2.set_title('(b) SPD')
ax2.set_xlabel('Longitude')
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig2, ax=ax2,orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# panel c - dust optical depth difference
fig3 = ax3.contourf(DDS.lon, DDS.lat, DDS.taudust_MY28_Ls270-DDS.taudust_MY24_Ls270,
                    levels=levels_dod_diff,cmap=py.cm.RdBu_r)
ax3.set_title('(c) Dust Optical Depth')
ax3.set_xlabel('Longitude')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
py.colorbar(fig3, ax=ax3,orientation='horizontal', ticks=[-5,-2.5,0,2.5,5])

py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300)

