# Figure S3 : Python script for base analysis and figures for Hartwick+2022
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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS3.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S3
# Enercon E33 (blue line), Jacobs 31-20 (orange line, scaled by a factor of 10),
# and NREL 5MW (green line, scaled by a factor of 1/10) compared with a
# theoretical power curve based on Eqn. 2 (purple curve). White bars (bottom)
# show annual, global average wind speed distribution. Filled bars (top) show the
# density-adjusted global average wind speed distribution used as inputs for
# turbine power curves.
# -------------------------------------------------------------------------------
ct = 3
print('Making Figure S3')

# Define Power Curves
bins = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13,13.5,14.,14.5,15,
        15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27.,27.5]
bins_5mw =  [0,1,2,3,4,5,6,7.,8.,9.,10.,11.,12.,13,14.,15,16,17,18,19,20,21,22,23,24,25,26,27.,28]
tmp = np.array(bins)
tmp_5mw = np.array(bins_5mw)

power1kwopt = np.zeros(len(bins))
power20kwopt = np.zeros(len(bins))
power333kwopt = np.zeros(len(bins))
power5mwopt = np.zeros(len(bins_5mw))

# general 1m2 turbine
S_1m = np.pi
cp_1m = 0.47
gen1m_powercurve = 1/2 * 0.0137 * tmp**3 * cp_1m * S_1m

# power for each turbine
thresh20kw = (0.,3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,11.,11.5, 12., 12.5,13.,13.5)
thresh20kw_top = (3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,11.,11.5, 12., 12.5,13.,13.5,14.)
power20kw = (0.,0.09,0.26,0.5,0.83,1.22,1.88,2.6,3.43,4.5,5.6,6.79,8.0,9.65,11.45,13.1,15.,17.2,18.7,19.8,20,20,20.)

power333kw = (0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
thresh333kw = (0.,2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13.)
thresh333kw_top = (2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 28.)

power5mw = (0.,40.5, 177.7, 403.9, 737.6, 1187.2, 1771.1, 2518.6, 3448.4, 4562.5, 5000.)
thresh5mw = (0.,3., 4., 5., 6., 7., 8., 9., 10., 11., 12.)
thresh5mw_top = (3.,4., 5., 6., 7., 8., 9., 10., 11., 25.)

#general E333 turbine
S_e333 = np.pi * (33.4/2)**2
cp_e333 = 0.3
genE333_powercurve = 1/2 * 1.125 * tmp**3 * cp_e333 * S_e333    # use 1.125 because we are comparing with adjusted wind speeds
#genE333_powercurve = 1/2 * 0.0137 * tmp**3 * cp_e333 * S_e333

# 20kW turbine
power20kwopt[6:28]=power20kw[1:]

# 333kW turbine
power333kwopt[4:27]=power333kw[:]
power333kwopt[27:]=335

# 5 MW turbine
power5mwopt[2:13]=power5mw[:]
power5mwopt=np.where(np.logical_and(tmp_5mw>=12, tmp_5mw<=25),5000,power5mwopt)

# Make Figure
fig, ax1 = py.subplots(1,1, figsize=(12,10))
fig.subplots_adjust(wspace=0.05)

fig1 = ax1.bar(bins, DDS_MY24.percent_time2_SY1*100, width=0.5,color='white', edgecolor='black')
fig1 = ax1.bar(bins,DDS_MY24.percent_time_SY1*100, width=0.5,hatch='/', color='grey', edgecolor='black',
        bottom=DDS_MY24.percent_time2_SY1*100)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

fig1 = ax2.plot(bins,power333kwopt, linewidth=3, label='Enercon E33, 330kW')
ax2.set_ylim(0,750)
ax2.plot(bins,genE333_powercurve*1e-3, linewidth=3, color='purple', label='30% Efficiency')
ax2.plot(bins,power20kwopt*10, linewidth=3, label='Jacobs 31-20, 20kW x 10')
ax2.plot(bins_5mw,power5mwopt*0.1, linewidth=3, label='NREL 5MW x 0.1')
ax2.set_ylabel('Turbine Power Output [kW]')
ax1.set_ylabel('Wind Speed Distribution [%]')
ax1.set_xlabel('Wind Speed [m/s]')
#ax2.legend()

angle=78
ax2.text(17.5,510, 'NREL 5MW x 0.1', color='tab:green')
ax2.text(17.5,345, 'Enercon E33', color='tab:blue')
ax2.text(12.5,210, 'Jacobs 31-20 x 10', color='tab:orange')
ax2.text(14.5,525,'30% Efficiency Turbine',rotation=angle, color='purple')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

