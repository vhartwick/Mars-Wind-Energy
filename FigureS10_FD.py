# Figure S10: Python script for base analysis and figures for Hartwick+2022
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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS10.nc'
DDS = xr.open_dataset(dataDIR, decode_times=False)

#-------------------------------------------------------------------------------
# FIGURE S10
# Enercon E33 Load duration curves (a) averaged globally and over the year and
# seasons (b) averaged annually at three locations of interest (c) at Protonilus
# Mensae (38N,48E) and (d) at Ismenius Lacus (29N,17E). The annual average lines
# in panels c and d reproduce the location annual average in panel b.The y-axis
# upper limit of panel a is 20% versus 80% in panels b-d.
# ------------------------------------------------------------------------------
ct = 10
print('Making Figure S10')

# 6. Plot
fig, axs = py.subplots(2,2,figsize=(12,5),sharex=True)
fig.subplots_adjust(wspace=0.1)
axs = axs.ravel()

# panel a global annual and seasonal average
fig = axs[0].plot(DDS.percent_time_gave_SY1,DDS.percent_max_gave, linewidth=2, label='Annual Average')

axs[0].plot(DDS.percent_time_ls0_SY1, DDS.percent_max_gave, linewidth=2, color='green',label='$\ L_{s}=0\degree$')
axs[0].plot(DDS.percent_time_ls90_SY1, DDS.percent_max_gave, linewidth=2, label='$\ L_{s}=90\degree$',color='orange')
axs[0].plot(DDS.percent_time_ls180_SY1, DDS.percent_max_gave, linewidth=2, label='$\ L_{s}=180\degree$',color='red')
axs[0].plot(DDS.percent_time_ls270_SY1, DDS.percent_max_gave, linewidth= 2, label='$\ L_{s}=270\degree$',color='purple')
axs[0].set_ylim([0,20])
axs[0].set_xlim([0.,100])
axs[0].grid(b=None,which='major', axis='both')
axs[0].set_title('(a) Global Average')
axs[0].set_ylabel('% Capacity Factor')
axs[0].legend()

# panel b - annual average at specific roi
fig = axs[1].plot(DDS.percent_time_gave_SY1,DDS.percent_max_gave, linewidth=2, label='Global Average')
axs[1].plot(DDS.percent_time_roi1_SY1,DDS.percent_max, 'red',linewidth=2, label='Deuteronilus Mensae')
axs[1].plot(DDS.percent_time_roi3_SY1, DDS.percent_max, 'purple',linewidth=2, label='Ismenius Lacus')
axs[1].plot(DDS.percent_time_roi4_SY1, DDS.percent_max, 'green', linewidth=2, label='Protonilus Mensae')
axs[1].set_ylim([0,80])
axs[1].grid(b=None,which='major', axis='both')
axs[1].set_title('(b) Annual Average')
axs[1].legend()

# panel c - Prontonilus Mensae
print('PM annual=', percent_time_roi4_SY1, percent_max)
print('PM Ls0=', ls0_PM, percent_max_PM)
print('PM Ls180=', ls180_PM, percent_max_PM)
print('PM Ls270=', ls270_PM, percent_max_PM)
fig = axs[2].plot(DDS.percent_time_roi4_SY1, DDS.percent_max, label='Annual Average')
axs[2].plot(DDS.ls0_PM, DDS.percent_max_PM, 'green', label='$\ L_{s}=0\degree$')
axs[2].plot(DDS.ls90_PM, DDS.percent_max_PM, 'orange', label='$\ L_{s}=90\degree$')
axs[2].plot(DDS.ls180_PM,DDS.percent_max_PM, 'red', label='$\ L_{s}=180\degree$')
axs[2].plot(DDS.ls270_PM, DDS.percent_max_PM, 'purple', label = '$\ L_{s}=270\degree$')
axs[2].set_ylim([0,80])
axs[2].grid(visible=None,which='major', axis='both')
axs[2].set_title('(c) Protonilus Mensae')
axs[2].set_xlabel('% of Time')
axs[2].set_ylabel('% Capacity Factor')
axs[2].legend()

# panel d - Ismenius Lacus
fig = axs[3].plot(DDS.percent_time_roi3_SY1, DDS.percent_max, label='Annual Average')
axs[3].plot(DDS.ls0_IL, DDS.percent_max_PM, 'green', label='$\ L_{s}=0\degree$')
axs[3].plot(DDS.ls90_IL, DDS.percent_max_PM, 'orange', label='$\ L_{s}=90\degree$')
axs[3].plot(DDS.ls180_IL,DDS.percent_max_PM, 'red', label='$\ L_{s}=180\degree$')
axs[3].plot(DDS.ls270_IL, DDS.percent_max_PM, 'purple', label = '$\ L_{s}=270\degree$')
axs[3].set_ylim([0,80])
axs[3].grid(visible=None,which='major', axis='both')
axs[3].set_title('(d) Ismenius Lacus')
axs[3].set_xlabel('% of Time')
axs[3].legend()

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)
