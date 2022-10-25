# Figure S8: Python script for base analysis and figures for Hartwick+2022
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
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/*.atmos_average.nc'
DDS_MY24 = xr.open_mfdataset(dataDIR,combine='by_coords',decode_times=False)

dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY28_highres/history/*.atmos_average.nc'
DDS_MY28 = xr.open_mfdataset(dataDIR,combine='by_coords',decode_times=False)


# -------------------------------------------------------------------------------
# FIGURE S8
# (a) Mars Year 24 and (b) MY 28 visible dust optical depth 
# -------------------------------------------------------------------------------
ct = 8
print('Making Figure S8')

# Plot
fig, ax = py.subplots(1,2, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)
levels = np.arange(0,3.6,0.1)

fig1 = ax[0].contourf(DDS_MY24.areo[:,0].squeeze()-360*3,DDS_MY24.lat,DDS_MY24.taudust_VIS.mean('lon').transpose(),levels=levels, cmap=py.cm.plasma)

ax[0].set_title('(a) MY24')
ax[0].set_ylabel('Latitude')
ax[0].set_xlabel('$\ L_{s}$')
ax[0].xaxis.set_major_locator(MultipleLocator(60))
ax[0].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig1,ax=ax[0])
cbar.set_label(r'[$\ \tau$]')

fig2 = ax[1].contourf(DDS_MY28.areo[:,0].squeeze()-360*3, DDS_MY24.lat, DDS_MY28.taudust_VIS.mean('lon').transpose(),levels=levels, cmap=py.cm.plasma)
ax[1].set_title('(b) MY28')
ax[1].set_xlabel('$\ L_{s}$')
ax[1].xaxis.set_major_locator(MultipleLocator(60))
ax[1].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig2,ax=ax[1])
cbar.set_label(r'[$\ \tau$]')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

#Send out Data
NEW_DF = xr.Dataset({'taudust_VIS_MY24':DDS_MY24.taudust_VIS,'areo_MY24':DDS_MY24.areo,
	'taudust_VIS_MY28':DDS_MY28.taudust_VIS,'areo_MY28':DDS_MY28.areo})
NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')

