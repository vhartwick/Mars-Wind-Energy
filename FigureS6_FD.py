# Figure S6: Python script for base analysis and figures for Hartwick+2022
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
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS6.nc'
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
# FIGURE S6
# Regions that produce (a,b) diurnal average and (c,d) nighttime average power
# levels greater than 24 kW (red squares), or greater than 14.2 kW (blue triangles)
# based on output from a single Enercon E33 at all sols/nights in the season
# (e.g. between Ls=80-100).
# -------------------------------------------------------------------------------
ct = 6
print('Making Figure S6')

# plot all locations with diurnal average power greater than some limit (ALL THE TIME)
# choose locations with energy (kW) always greater than the minimum energy requirements (24)
# or the average greater than the mim energy requirement

# 3. Plot
# SEASONAL DIURNAL AVERAGE
fig, [(ax1,ax2),(ax3,ax4)] = py.subplots(2,2,sharey=True, sharex=True,figsize=(15,8))
py.subplots_adjust(wspace=0.05)

# panel a - Ls90
fig1 = ax1.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf/1000, 20,cmap=py.cm.Greys_r)
clb=py.colorbar(fig1,ax=ax1)
clb.set_label('[km]')
ax1.set_title('(a) $\ L_{s}=90\degree$')
ax1.set_ylabel('Latitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
# add landing sites
ax1.plot(longitude, latitude, 'o',markerfacecolor='green',label='ROIs')
#ax1.plot(DDS_MY24.lon_sav2_ls90, DDS_MY24.lat_sav2_ls90, 'ok', markersize=5,fillstyle='none', label='>2.2kW')
ax1.plot(DDS_MY24.lon_sav15_ls90, DDS_MY24.lat_sav15_ls90, 'vk', markersize=10, markerfacecolor='blue',label='>14.2kW')
ax1.plot(DDS_MY24.lon_sav24_ls90, DDS_MY24.lat_sav24_ls90, 'sk', markersize=10, markerfacecolor='red', label='>24kW')

# panel b - Ls270
fig2 = ax2.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf/1000, 20,cmap=py.cm.Greys_r)
clb=py.colorbar(fig2,ax=ax2)
clb.set_label('[km]')
ax2.set_title('(b) $\ L_{s}=270\degree$')
# add landing sites
ax2.plot(longitude, latitude, 'o',markerfacecolor='green',label='ROIs')
#ax2.plot(DDS_MY24.lon_sav2_ls270, DDS_MY24.lat_sav2_ls270, 'ok', markersize=5,fillstyle='none', label='>2.2kW')
ax2.plot(DDS_MY24.lon_sav15_ls270, DDS_MY24.lat_sav15_ls270, 'vk', markersize=10, markerfacecolor='blue',label='>14.2kW')
ax2.plot(DDS_MY24.lon_sav24_ls270, DDS_MY24.lat_sav24_ls270, 'sk', markersize=10, markerfacecolor='red', label='>24kW')

# panel c - nighttime Ls90
fig3 = ax3.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf/1000, 20,cmap=py.cm.Greys_r)
clb=py.colorbar(fig3,ax=ax3)
clb.set_label('[km]')
ax3.set_title('(c) Night, $\ L_{s}=90\degree$')
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')
# add landing sites
ax3.plot(longitude, latitude,'o',markerfacecolor='green',label='ROIs')
#ax3.plot(DDS_MY24.lon_sav2_ls90pm, DDS_MY24.lat_sav2_ls90pm, 'ok', markersize=5,fillstyle='none', label='>2.2kW')
ax3.plot(DDS_MY24.lon_sav15_ls90pm, DDS_MY24.lat_sav15_ls90pm, 'vk', markersize=10, markerfacecolor='blue',label='>14.2kW')
ax3.plot(DDS_MY24.lon_sav24_ls90pm, DDS_MY24.lat_sav24_ls90pm, 'sk', markersize=10, markerfacecolor='red', label='>24kW')

# panel d - night, Ls270
fig4 = ax4.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS_MY24.zsurf/1000., 20,cmap=py.cm.Greys_r)
clb=py.colorbar(fig4,ax=ax4)
clb.set_label('[km]')
ax4.set_title('(d) Night, $\ L_{s}=270\degree$')
ax4.set_xlabel('Longitude')

# add landing sites
ax4.plot(longitude, latitude,'o',markerfacecolor='green',label='ROIs')
#ax4.plot(DDS_MY24.lon_sav2_ls270pm, DDS_MY24.lat_sav2_ls270pm, 'ok', markersize=5,fillstyle='none', label='>2.2kW')
ax4.plot(DDS_MY24.lon_sav15_ls270pm, DDS_MY24.lat_sav15_ls270pm, 'vk', markersize=10, markerfacecolor='blue',label='>14.2kW')
ax4.plot(DDS_MY24.lon_sav24_ls270pm, DDS_MY24.lat_sav24_ls270pm, 'sk', markersize=10, markerfacecolor='red', label='>24kW')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

