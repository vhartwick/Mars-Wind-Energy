# Figure S11: Python script for base analysis and figures for Hartwick+2022
# Author: Victoria Hartwick
# Last Modified: 6/29/22

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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS11.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

#define prospective landing sites of interest
longitude = [2, 22.117202268431015, 22.797731568998103, 17.0132325141777, 47.97731568998114, 55.4631379962193,
             76.89981096408323, 76.89981096408323, 77.58034026465029, 96.63516068052931, 94.25330812854443,
             102.07939508506618,104.12098298676747, 126.23818525519849,137.80718336483935, 172.17391304347828,
             162.30623818525518, 158.22306238185254, 176.93761814744803,174.21550094517963,175.9168241965974,
             188.83636363636364,201.92727272727274, 194.07272727272726,199.96363636363637,190.8,220.25454545454545,
             267.3818181818182,270.9818181818182, 287.3454545454546, 292.90909090909093,289.9636363636364,
             300.1090909090909,311.5636363636364, 320.0727272727273, 313.8545454545455, 325.30909090909097,
             310.9090909090909,323.0181818181818, 341.0181818181818, 337.41818181818184, 338.72727272727275,
             351.81818181818187, 353.4545454545455,355.41818181818184,350.5090909090909] 

latitude = [0, 36.702439024390266, 34.75609756097563, 29.195121951219527, 37.81463414634148, -11.121951219512184,
            18.62926829268295, 14.18048780487807, -23.07804878048779, -27.526829268292673, -31.41951219512194,
            -34.756097560975604, -35.868292682926835, 16.68292682926831, -3.6146341463414444, 34.75609756097563,
            8.063414634146355, -0.834146341463395, -10.009756097560967, -10.843902439024376, -12.234146341463415,
            42.5831381733021, 34.840749414519905, -24.962529274004694, -36.17564402810305, -46.32084309133489,
            35.64168618266979, -5.473067915690862, 38.04449648711944, 21.224824355971897, -9.477751756440284,
            -10.011709601873534, -25.229508196721312, 19.088992974238877, 16.686182669789233, 9.744730679156902,
            -3.3372365339578494, -16.4192037470726, -31.370023419203747, 20.69086651053864, 18.555035128805613,
            2.002341920374704, -0.13348946135830886, -1.46838407494144, -2.269320843091336, -37.777517564402814]

# locate indeces of ROI lon/lat
lon_idx, lat_idx = [], []
tmp = np.array(longitude)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lon-tmp[i])).argmin()
    lon_idx += [idx]

tmp = np.array(latitude)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_MY24.lat-tmp[i])).argmin()
    lat_idx += [idx]

loc_name = ['Meridian Planum', 'Deuteronilus Mensae', 'DM2', 'Ismenius Lacus', 'Protonilus Mensae', 'Huygens Crater',
            'Nili Fossae', 'Jezero Crater', 'Hadriacus Palus', 'Ausonia Cavus', 'Hellas Rim', 'Mesopotamia',
            'Eastern Hellas', 'Hebrus Valles', 'Gale Crater', 'Phlegra Dorsa', 'Cerberus',
            'Hills Zephyria Planum', 'Apollinaris Sulci', 'AS2', 'Guzev Crater', 'Amazonis Planitia', 'Erberus Montes',
            'Columbus Crater', 'Newton Crater', 'Copernicus Crater', 'Acheron Fossae', 'Noctis Landing', 'Tempe Terra',
           'Kasei Valles', 'Coprates Chasma', 'Melas Casma', 'Southern Nectaris Fossae', 'Chryse/Viking',
            'Valles Marineris Mouth', 'Hypanis Valles', 'Eastern Valles Marineris', 'Equatorial Valles Marineris',
            'Hale Crater', 'Mawrth Valles',
            'McLaughlin Crater', 'Aram Chaos', 'Firsoff Crater', 'Sinus Meridiani', 'Southern Meridiani', 'Noachis Terra']

# -------------------------------------------------------------------------------
# FIGURE S11
# (a,d) Annual average power production for a Jacobs 31-20 and NREL 5MW wind
# turbine, (b,e) the additional percent of the Mars year that power generation
# exceeds 24kW for a solar array (Cp=0.2, S=2500m2) and a single turbine, and
# (c,f) regions that produce diurnal average power greater than 0.2kW (blue
# triangles) and 0.5 kW (red squares) for the Jacobs 31-20 turbine and greater
# than 14.2 kW (blue triangles) or 24kW (red squares) for the NREL 5MW turbine at
# all times throughout the Mars year.
# -------------------------------------------------------------------------------
ct = 11
print('Making Figure S11')

# 3. Plot
energylim=[0,24,1000]  # energy requirement per sol (Rucker 2015 https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20150000526.pdf)
energylim_JAC = [0, 1, 1000]
slevels = np.arange(0,1600,100)
smin = np.arange(0,125,25)
smin2 = np.arange(10000,100000,20000)

fig, ax = py.subplots(2,3,figsize=(12,8),sharex=True, sharey=True)
fig.subplots_adjust(wspace=0.13)
ax = ax.ravel()

# Jacobs
fig = ax[0].contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.aep_Jac[:,:], 20, cmap=py.cm.viridis)
ax[0].set_ylabel('Latitude')
ax[0].set_title('(a) AEP')
ax[0].xaxis.set_major_locator(MultipleLocator(60))
ax[0].yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig,ax=ax[0])
clb.set_label('[kWh]')

fig = ax[1].contourf(DDS_MY24.lon, DDS_MY24.lat,DDS_MY24.percent_Jac*100-DDS_PCT.percent_solar2500*100,cmap=py.cm.magma)
clb=py.colorbar(fig,ax= ax[1])
clb.set_label('[% time]')
py.gcf().text(0.05, 0.65, 'Jacobs 31-20', fontsize=12, rotation=90)
ax[1].xaxis.set_ticklabels([])
ax[1].set_title('(b) % Increase')

fig = ax[2].contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf/1000., 10,cmap=py.cm.Greys_r)
clb=py.colorbar(fig,ax= ax[2])
clb.set_label('[km]')
ax[2].xaxis.set_major_locator(MultipleLocator(60))
ax[2].yaxis.set_major_locator(MultipleLocator(30))
ax[2].plot(longitude, latitude, 'o',markerfacecolor='green',label='ROIs')
ax[2].plot(DDS_MY24.lon_sav15_JAC, DDS_MY24.lat_sav15_JAC, 'vk', markersize=8, markerfacecolor='blue', label='>0.2kW')
ax[2].plot(DDS_MY24.lon_sav24_JAC, DDS_MY24.lat_sav24_JAC, 'sk', markersize=8, markerfacecolor='red', label='>0.5kW')
#ax[5].legend()
ax[2].set_title('(c) Available Power')

# NREL
fig = ax[3].contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.aep_NREL[:,:], 20, cmap=py.cm.viridis)
ax[3].set_ylabel('Latitude')
ax[3].set_xlabel('Longitude')
ax[3].set_title('(d) AEP')
ax[3].xaxis.set_major_locator(MultipleLocator(60))
ax[3].yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig,ax=ax[3])
clb.set_label('[GWh]')

fig = ax[4].contourf(DDS_MY24.lon, DDS_MY24.lat,DDS_MY24.percent_NREL*100-DDS_PCT.percent_solar2500*100,cmap=py.cm.magma)
ax[4].set_xlabel('Longitude')
clb=py.colorbar(fig,ax= ax[4])
clb.set_label('[% time]')
ax[4].set_title('(e) % Increase')
py.gcf().text(0.05, 0.23, 'NREL 5MW', fontsize=12, rotation=90)

fig = ax[5].contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf/1000., 10,cmap=py.cm.Greys_r)
ax[5].set_xlabel('Longitude')
clb=py.colorbar(fig,ax= ax[5])
ax[5].set_title('(f) Available Power')
clb.set_label('[km]')
ax[5].xaxis.set_major_locator(MultipleLocator(60))
ax[5].yaxis.set_major_locator(MultipleLocator(30))
ax[5].plot(longitude, latitude, 'o',markerfacecolor='green',label='ROIs')
ax[5].plot(DDS_MY24.lon_sav15_NREL, DDS_MY24.lat_sav15_NREL, 'vk', markersize=2, markeredgecolor='blue',markerfacecolor='blue', label='>15kW')
ax[5].plot(DDS_MY24.lon_sav24_NREL, DDS_MY24.lat_sav24_NREL, 'sk', markersize=2,markeredgecolor='red', markerfacecolor='red', label='>24kW')
#ax[8].legend()

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

