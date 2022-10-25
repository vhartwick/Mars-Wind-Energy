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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS9.nc'
DDS = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')

ct = 9 
print('Making Figure S9')

# Plot and Save Figure
fig, ax = py.subplots()
ax.set_title('Wind Speeds at Insight Landing Site ($\ 4.5\degree N,135\degree E$)')
ax.set_ylabel('[m/s]')
ax.set_xlabel('$\ L_{s}$')
fig = py.plot(DDS.areo[:,0,0].squeeze()-360*3, DDS.ws_twins.mean('time_of_day_24'),"s")
py.plot(DDS.areo_i,DDS.ws_i,"^")
ax.set_xlim([0,360])
ax.xaxis.set_major_locator(MultipleLocator(60))

ax.text(20,5, 'Observed  Winds', color='tab:orange')
ax.text(215,2.75, 'Simulated Winds', color='tab:blue')
py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)
