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
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_pow.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

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

# 1. Define Seasons
yrs = [0,1,2,3]
areo = DDS_MY24.areo[:,0]
lsrange, lsrange1, lsrange2, lsrange3 = [], [], [], []
for i in yrs:
    loc = np.where(np.logical_and(areo >= 360*i + 45, areo <= 360*i + 135.))
    loc1 = np.where(np.logical_and(areo >= 360*i + 135, areo <= 360*i + 225.))
    loc2 = np.where(np.logical_and(areo >= 360*i + 225, areo <= 360*i + 315.))
    loc3 = np.where(np.logical_and(areo >= 360*i + 315, areo <= 360*(i+1) + 45.))

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

# 2. Calculate Load Duration Curves for Global Year & Season
# calculate based on power (preferred)
p333 = (0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
p333top = (3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
power_range =(0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)

# calculate global mean 
zonal_ave_SY1 = np.mean(DDS_MY24.pow_diurn_E33, axis=3)

# make sure to weight by area
latr = np.deg2rad(DDS_MY24.lat)
weights = np.cos(latr)
global_ave_SY1 = np.average(zonal_ave_SY1, axis=2, weights=weights)

# select seasons
global_ave_ls90_SY1 = global_ave_SY1[ls90,:]
global_ave_ls180_SY1 = global_ave_SY1[ls180,:]
global_ave_ls270_SY1 = global_ave_SY1[ls270,:]
global_ave_ls0_SY1 = global_ave_SY1[ls0,:]

# what we really want is the amount of time capacity factor is x or greater so...
count_SY1, countls90_SY1, countls180_SY1, countls270_SY1, countls0_SY1 = [], [], [], [], []

for x in power_range:   # loop through

    # SY1
    count_SY1 += [np.count_nonzero(np.where(global_ave_SY1[:,:] >= x))]
    countls90_SY1 += [np.count_nonzero(np.where(global_ave_ls90_SY1[:,:] >= x))]
    countls180_SY1 += [np.count_nonzero(np.where(global_ave_ls180_SY1[:,:] >= x))]
    countls270_SY1 += [np.count_nonzero(np.where(global_ave_ls270_SY1[:,:] >= x))]
    countls0_SY1 += [np.count_nonzero(np.where(global_ave_ls0_SY1[:,:] >= x))]

# calculate percent of time operating in each interval
percent_time_gave_SY1 = np.divide(count_SY1,count_SY1[0])*100
percent_time_ls90_SY1 = np.divide(countls90_SY1,countls90_SY1[0])*100
percent_time_ls180_SY1 = np.divide(countls180_SY1,countls180_SY1[0])*100
percent_time_ls270_SY1 = np.divide(countls270_SY1,countls270_SY1[0])*100
percent_time_ls0_SY1 = np.divide(countls0_SY1,countls0_SY1[0])*100

# calculate capacity factor (power interval / total)
percent_max_gave = np.divide(power_range, max(power_range))*100
print(percent_max_gave)
print(percent_time_gave_SY1)

# 3. Now calculate Load Duration Curves (Annual Average) at specific ROI
lon_roi = [22.117202268431015, 22.797731568998103, 17.0132325141777,47.97731568998114]
lat_roi = [36.702439024390266, 34.75609756097563, 29.195121951219527,37.81463414634148]
roi_name = ['Deuteronilus Mensae', ' ', 'Ismenius Lacus', 'Protonilus Mensae']
color = ['red','purple','green']

# locate indeces of ROI lon/lat
lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lon-tmp[i])).argmin()
    lonlim_idx += [idx]

tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lat-tmp[i])).argmin()
    latlim_idx += [idx]

local_power_SY1 = DDS_MY24.pow_diurn_E33

# count number of times turbine is operating at or above power levels
count_roi1_SY1,count_roi2_SY1,count_roi3_SY1,count_roi4_SY1 = [],[],[],[]

for x in power_range:   # loop through
    tmp_roi1_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_roi1_SY1 += [tmp_roi1_SY1]


    tmp_roi2_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[1],lonlim_idx[1]] >= x))
    count_roi2_SY1 +=[tmp_roi2_SY1]


    tmp_roi3_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[2],lonlim_idx[2]] >= x))
    count_roi3_SY1 +=[tmp_roi3_SY1]

    tmp_roi4_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[3],lonlim_idx[3]] >= x))
    count_roi4_SY1 +=[tmp_roi4_SY1]

# calculate percent of time turbine is operating at each power interval
percent_time_roi1_SY1 = np.divide(count_roi1_SY1, count_roi1_SY1[0])*100
percent_time_roi2_SY1 = np.divide(count_roi2_SY1, count_roi2_SY1[0])*100
percent_time_roi3_SY1 = np.divide(count_roi3_SY1, count_roi3_SY1[0])*100
percent_time_roi4_SY1 = np.divide(count_roi4_SY1, count_roi4_SY1[0])*100

# calculate capacity factor (power interval / total)
percent_max = np.divide(power_range, max(power_range))*100

# 4. Deep Dive @ Prontonilus Mensae - Look at One Location for all Seasons
lon_roi = [47.97731568998114]
lat_roi = [37.81463414634148]
roi_name = ['Protonilus Mensae']

# locate indeces of ROI lon/lat
lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lon-tmp[i])).argmin()
    lonlim_idx += [idx]

tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lat-tmp[i])).argmin()
    latlim_idx += [idx]

local_power_SY1 = DDS_MY24.pow_diurn_E33

# this plots the power versus time at that location
# count number of times turbine is operating at or above power levels

count_ls90_SY1,count_ls180_SY1,count_ls270_SY1,count_ls0_SY1 = [],[],[],[]
for x in power_range:   # loop through
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls90,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls90_SY1 += [tmp_SY1]

    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls180,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls180_SY1 +=[tmp_SY1]

    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls270,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls270_SY1 +=[tmp_SY1]


    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls0,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls0_SY1 +=[tmp_SY1]

# calculate percent of time turbine is operating at each power interval
ls90_PM = np.divide(count_ls90_SY1, count_ls90_SY1[0])*100
ls180_PM = np.divide(count_ls180_SY1, count_ls180_SY1[0])*100
ls270_PM = np.divide(count_ls270_SY1, count_ls270_SY1[0])*100
ls0_PM = np.divide(count_ls0_SY1, count_ls0_SY1[0])*100

# calculate capacity factor (power interval / total)
percent_max_PM = np.divide(power_range, max(power_range))*100

#5. Finally, Ismenius Lacus
# ismenius lacus quadrangle 0-60E, 30-65 N
lon_roi = [17.0132325141777]
lat_roi = [29.195121951219527]
roi_name = ['Ismenius Lacus']

# locate indeces of ROI lon/lat
lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lon-tmp[i])).argmin()
    lonlim_idx += [idx]

tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lat-tmp[i])).argmin()
    latlim_idx += [idx]

local_power_SY1 = DDS_MY24.pow_diurn_E33

# this plots the power versus time at that location
# count number of times turbine is operating at or above power levels

count_ls90_SY1,count_ls180_SY1,count_ls270_SY1,count_ls0_SY1 = [],[],[],[]
for x in power_range:   # loop through
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls90,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls90_SY1 += [tmp_SY1]

    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls180,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls180_SY1 +=[tmp_SY1]

    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls270,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls270_SY1 +=[tmp_SY1]


    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls0,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls0_SY1 +=[tmp_SY1]

# calculate percent of time turbine is operating at each power interval
ls90_IL = np.divide(count_ls90_SY1, count_ls90_SY1[0])*100
ls180_IL = np.divide(count_ls180_SY1, count_ls180_SY1[0])*100
ls270_IL = np.divide(count_ls270_SY1, count_ls270_SY1[0])*100
ls0_IL = np.divide(count_ls0_SY1, count_ls0_SY1[0])*100

# calculate capacity factor (power interval / total)
percent_max_IL = np.divide(power_range, max(power_range))*100

# 6. Plot
fig, axs = py.subplots(2,2,figsize=(12,5),sharex=True)
fig.subplots_adjust(wspace=0.1)
axs = axs.ravel()

# panel a global annual and seasonal average
fig = axs[0].plot(percent_time_gave_SY1,percent_max_gave, linewidth=2, label='Annual Average')

axs[0].plot(percent_time_ls0_SY1, percent_max_gave, linewidth=2, color='green',label='$\ L_{s}=0\degree$')
axs[0].plot(percent_time_ls90_SY1, percent_max_gave, linewidth=2, label='$\ L_{s}=90\degree$',color='orange')
axs[0].plot(percent_time_ls180_SY1, percent_max_gave, linewidth=2, label='$\ L_{s}=180\degree$',color='red')
axs[0].plot(percent_time_ls270_SY1, percent_max_gave, linewidth= 2, label='$\ L_{s}=270\degree$',color='purple')
axs[0].set_ylim([0,20])
axs[0].set_xlim([0.,100])
axs[0].grid(b=None,which='major', axis='both')
axs[0].set_title('(a) Global Average')
axs[0].set_ylabel('% Capacity Factor')
axs[0].legend()

# panel b - annual average at specific roi
fig = axs[1].plot(percent_time_gave_SY1,percent_max_gave, linewidth=2, label='Global Average')
axs[1].plot(percent_time_roi1_SY1,percent_max, 'red',linewidth=2, label='Deuteronilus Mensae')
axs[1].plot(percent_time_roi3_SY1, percent_max, 'purple',linewidth=2, label='Ismenius Lacus')
axs[1].plot(percent_time_roi4_SY1, percent_max, 'green', linewidth=2, label='Protonilus Mensae')
axs[1].set_ylim([0,80])
axs[1].grid(b=None,which='major', axis='both')
axs[1].set_title('(b) Annual Average')
axs[1].legend()

# panel c - Prontonilus Mensae
print('PM annual=', percent_time_roi4_SY1, percent_max)
print('PM Ls0=', ls0_PM, percent_max_PM)
print('PM Ls180=', ls180_PM, percent_max_PM)
print('PM Ls270=', ls270_PM, percent_max_PM)
fig = axs[2].plot(percent_time_roi4_SY1, percent_max, label='Annual Average')
axs[2].plot(ls0_PM, percent_max_PM, 'green', label='$\ L_{s}=0\degree$')
axs[2].plot(ls90_PM, percent_max_PM, 'orange', label='$\ L_{s}=90\degree$')
axs[2].plot(ls180_PM,percent_max_PM, 'red', label='$\ L_{s}=180\degree$')
axs[2].plot(ls270_PM, percent_max_PM, 'purple', label = '$\ L_{s}=270\degree$')
axs[2].set_ylim([0,80])
axs[2].grid(visible=None,which='major', axis='both')
axs[2].set_title('(c) Protonilus Mensae')
axs[2].set_xlabel('% of Time')
axs[2].set_ylabel('% Capacity Factor')
axs[2].legend()

# panel d - Ismenius Lacus
fig = axs[3].plot(percent_time_roi3_SY1, percent_max, label='Annual Average')
axs[3].plot(ls0_IL, percent_max_PM, 'green', label='$\ L_{s}=0\degree$')
axs[3].plot(ls90_IL, percent_max_PM, 'orange', label='$\ L_{s}=90\degree$')
axs[3].plot(ls180_IL,percent_max_PM, 'red', label='$\ L_{s}=180\degree$')
axs[3].plot(ls270_IL, percent_max_PM, 'purple', label = '$\ L_{s}=270\degree$')
axs[3].set_ylim([0,80])
axs[3].grid(visible=None,which='major', axis='both')
axs[3].set_title('(d) Ismenius Lacus')
axs[3].set_xlabel('% of Time')
axs[3].legend()

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

#Send out Data
NEW_DF = xr.Dataset({'percent_max_gave':percent_max_gave,'percent_max':percent_max,'percent_max_PM':percent_max_PM,
	'percent_time_gave_SY1':percent_time_gave_SY1,
	'percent_time_ls0_SY1':percent_time_ls0_SY1,'percent_time_ls90_SY1':percent_time_ls90_SY1,
	'percent_time_ls180_SY1':percent_time_ls180_SY1,'percent_time_ls270_SY1':percent_time_ls270_SY1,
	'percent_time_roi1_SY1':percent_time_roi1_SY1,'percent_time_roi2_SY1':percent_time_roi2_SY1,'percent_time_roi1_SY3':percent_time_roi3_SY1,
	'percent_time_roi4_SY1':percent_time_roi4_SY1,'ls0_PM':ls0_PM,'ls90_PM':ls90_PM,'ls180_PM':ls180_PM,'ls270_PM':ls270_PM,
	'ls0_IL':ls0_IL,'ls90_IL':ls90_IL,'ls180_IL':ls180_IL,'ls270_IL':ls270_IL})
NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')
