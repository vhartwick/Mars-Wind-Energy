# python routine to plot environment at DM2 (35N,23E)
# author: victoria l. hartwick
# last modified: 6/9/22


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
from windrose import WindroseAxes
import pandas as pd


# Get Path
PATH = os.getcwd()

# Import Base Datasets
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/*.atmos_diurn_T_zagl_min.nc'
DDS_MY24 = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')

ct = 9 
print('Making Figure S9')
u = DDS_MY24.ucomp.sel(zagl=5).sel(lat=4.5,method='nearest').sel(lon=136,method='nearest')
v = DDS_MY24.vcomp.sel(zagl=5).sel(lat=4.5,method='nearest').sel(lon=136,method='nearest')
ws = (u**2 + v**2)**0.5

# Interpolate to 1.5m
zo = 0.01 # roughness length
u_twins = u * log(1.6/zo)/log(5/zo)
v_twins = v * log(1.6/zo)/log(5/zo)
ws_twins = (u_twins**2 + v_twins**2)**0.5

# Import Insight WS Data
dataDIR = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/Insight_savefile.nc'
DFI = xr.open_dataset(dataDIR,decode_times=False)

ws_i = np.array(DFI.ws.mean('tod'))
ws_i[31:34]=np.nan
ws_i[100:104]=np.nan
ws_i[120:124]=np.nan
print(DDS_MY24.areo)
print(ws_i)
#ws_i[ws_i==0]=np.nan
areo_i = DFI.insight_areo

# Plot and Save Figure
fig, ax = py.subplots()
ax.set_title('Wind Speeds at Insight Landing Site ($\ 4.5\degree N,135\degree E$)')
ax.set_ylabel('[m/s]')
ax.set_xlabel('$\ L_{s}$')
fig = py.plot(DDS_MY24.areo[:,0,0].squeeze()-360*3, ws_twins.mean('time_of_day_24'),"s")
py.plot(areo_i,ws_i,"^")
ax.set_xlim([0,360])
ax.xaxis.set_major_locator(MultipleLocator(60))

ax.text(20,5, 'Observed  Winds', color='tab:orange')
ax.text(215,2.75, 'Simulated Winds', color='tab:blue')
py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

#Send out Data
var1 = xr.DataArray(ws_twins, coords={'time':DDS_MY24.time, 'time_of_day_24':DDS_MY24.time_of_day_24},dims=['time','time_of_day_24'])
var1.name = 'ws_twins'

var2 = xr.DataArray(DDS_MY24.areo[:,0,0].squeeze(), coords={'time':DDS_MY24.time},dims=['time'])
var2.name = 'areo'

NEW_DF = xr.merge([var1,var2])
NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')
