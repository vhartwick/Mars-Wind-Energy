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
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY34/AAV_300W_power_MY34_LOWRES.nc'
DDS = xr.open_dataset(dataDIR, decode_times=False)

print(DDS)
# locate indeces of ROI lon/lat

lon_roi = [2.28]
lat_roi = [-5.23]

lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS.lon-tmp[i])).argmin()
    lonlim_idx += [idx]
        
tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS.lat-tmp[i])).argmin()
    latlim_idx += [idx]
    

# convert kW to W by multiplying by 1000
dust = DDS.taudust_VIS[77,:, 22, 0]*10
ls190_exact = DDS.pow_diurn_SY34[77,:,22,0] * 1000

# -------------------------------------------------------------------------------
# FIGURE S9
# -------------------------------------------------------------------------------
ct = 9
print('Making Figure S9')

py.xlabel('Time of Day (hours)')
py.ylabel('Wind Power Output (watts)')
py.axhline(1, color='red', linestyle='--', label='Energy Requirment')
py.plot(DDS.time_of_day_24, ls190_exact, color='purple')
py.plot(DDS.time_of_day_24, dust, color='orange', label='Dust Visibility')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)


