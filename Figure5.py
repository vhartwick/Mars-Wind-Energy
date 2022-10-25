# Figure 5 : Python script for base analysis and figures for Hartwick+2022
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

dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_savefile.nc'
DDS = xr.open_dataset(dataDIR, decode_times=False)

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
# FIGURE 5
# Regions that produce 5-sol diurnal average power levels continuously through a
# year greater than 24 kW (red squares), greater than 14 kW (blue triangles) or
# greater then 2.2kW (open circles) based on output from a single Enercon E33
# turbine. Solid black circles show proposed landing sites from the NASA HLS2
# Workshop. Contour lines show regions where the annual mean wind power production
# is greater than 14.2 and 24kW. Filled contours show the topography.
# -------------------------------------------------------------------------------
ct = 5
print('Making Figure 5')

# plot all locations with diurnal average power greater than some limit (ALL THE TIME)
# choose locations with energy (kW) always greater than the minimum energy requirements (24)
# or the average greater than the mim energy requirement

# diurnal mean
tmp_SY1 = np.mean(DDS_MY24.pow_diurn_E33,axis=1)

lat_sav24, lon_sav24 = [], []
lat_sav15, lon_sav15 = [], []
lat_sav2, lon_sav2 = [], []

for i in range(0,len(DDS_MY24.lat)):
    for j in range(0,len(DDS_MY24.lon)):
        tmp = tmp_SY1[:,i,j]
        if np.any(tmp < 24):
            continue
        else:
            lat_sav24+= [DDS_MY24.lat[i]]
            lon_sav24 +=[DDS_MY24.lon[j]]

for i in range(0,len(DDS_MY24.lat)):
    for j in range(0,len(DDS_MY24.lon)):
        tmp = tmp_SY1[:,i,j]
        if np.any(tmp < 14.2):
            continue
        else:
            lat_sav15 += [DDS_MY24.lat[i]]
            lon_sav15 +=[DDS_MY24.lon[j]]

for i in range(0,len(DDS_MY24.lat)):
    for j in range(0,len(DDS_MY24.lon)):
        tmp = tmp_SY1[:,i,j]
        if np.any(tmp < 2.2):
            continue
        else:
            lat_sav2 += [DDS_MY24.lat[i]]
            lon_sav2 +=[DDS_MY24.lon[j]]

energy_per_sol = np.sum(DDS_MY24.pow_diurn_E33,axis=1).mean('time')
tmp2 = DDS.wpd.mean(('time','time_of_day_24')).sel(zagl=50)
aep = np.sum(DDS_MY24.pow_diurn_E33,axis=0)
aep = np.sum(aep,axis=0)*5
for i in range(0,len(np.array(lat_sav24))):
	print(np.array(lat_sav24)[i], np.array(lon_sav24)[i],np.array(tmp2.sel(lat=np.array(lat_sav24)[i],lon=np.array(lon_sav24)[i],method='nearest')),np.array(energy_per_sol.sel(lat=np.array(lat_sav24)[i],lon=np.array(lon_sav24)[i],method='nearest')),np.array(aep.sel(lat=np.array(lat_sav24)[i],lon=np.array(lon_sav24)[i],method='nearest')))

# highly restrictive so also plot annual average power > 24
tmp = DDS_MY24.pow_diurn_E33.mean('time')
annual_ave_pow_E33 = tmp.mean('time_of_day_24')

# plot
level_lim = [0,14.2,24,1000]
colors =['tab:blue','tab:orange']
fig, ax4 = py.subplots(1,1,sharey=True, figsize=(15,8))

fig4 = ax4.contourf(DDS_MY24.lon, DDS_MY24.lat, DDS.zsurf/1000, 20,cmap=py.cm.Greys)
clb=py.colorbar(fig4,ax=ax4)
clb.set_label('[km]')
fig = ax4.contour(DDS_MY24.lon, DDS_MY24.lat, annual_ave_pow_E33, colors=colors, levels=level_lim)
py.clabel(fig, fontsize=9,inline=1, fmt='%1.0f')

ax4.plot(lon_sav2, lat_sav2, 'ok', markersize=10,fillstyle='none', label='>2.2kW')
ax4.plot(lon_sav15, lat_sav15, 'vk', markersize=15, markerfacecolor='blue',label='>14.2kW')
ax4.plot(lon_sav24, lat_sav24, 'sk', markersize=15, markerfacecolor='red', label='>24kW')
ax4.set_xlabel('Longitude')
ax4.set_ylabel('Latitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))

# add landing sites
ax4.plot(longitude, latitude, 'o', color='#02d8e9',markersize=12, markeredgecolor='k',markeredgewidth=1.5,label='ROIs')

ax4.legend(ncol=4, loc='lower center')
py.savefig(f'{PATH}/WindEnergy_HighRes_Fig{ct}.eps',dpi=300)

# Send out Data
var1 = xr.DataArray(annual_ave_pow_E33, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var1.name = 'annual_ave_pow_E33'

var2 = xr.DataArray(DDS.zsurf, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var2.name = 'zsurf'

NEW_DF = xr.merge([var1,var2])
NEW_DF.to_netcdf(f'{PATH}/Data/Figure{ct}.nc')

