# Python Routine to perform base analysis for Hartwick+2022 Wind Energy
# Author: Victoria Hartwick
# Last Updated: March 2, 2022

# Import Scientific Modules
import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import Dataset, MFDataset
from numpy import sqrt, exp, max, mean, min, log, log10,int
from itertools import product
import xarray as xr

# ========================================================
# Analysis of MY24
# ========================================================
# Open Data Files
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/02094.fixed.nc'
DS_fixed = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/*.atmos_diurn_T_zagl_min.nc'
DDS = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')

# 1. Calculate Wind Power Density
# calculate wind speed (again using diurnal file)
vel_diurn = sqrt(DDS.ucomp**2 + DDS.vcomp**2)

# calculate wind profile denisty [W/m2]
wpd = 0.5 * DDS.rho * vel_diurn**3

# calculate adjusted wind speed
utmp = DDS.ucomp * (DDS.rho/1.125)**(1./3.)
vtmp= DDS.vcomp * (DDS.rho/1.125)**(1./3.)
veladj= sqrt(utmp**2 + vtmp**2)

# Output as New NC File

# use atmospheric average file for short wave solar surface energy flux [W/m2]
swflx = DDS.swflx
rho = DDS.rho

var1 = xr.DataArray(swflx, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var1.name = 'swflx'

var2 = xr.DataArray(wpd, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24, 'zagl': DDS.zagl},dims=['time','time_of_day_24','zagl','lat','lon'])
var2.name = 'wpd'

var3 = xr.DataArray(veladj, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24, 'zagl': DDS.zagl},dims=['time','time_of_day_24','zagl','lat','lon'])
var3.name = 'veladj'

var4 = xr.DataArray(vel_diurn, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24, 'zagl': DDS.zagl},dims=['time', 'time_of_day_24','zagl','lat','lon'])
var4.name = 'ws'

var5 = xr.DataArray(DDS.areo, coords={'time': DDS.time,'time_of_day_24': DDS.time_of_day_24, 'scalar_axis': DDS.scalar_axis},dims=['time', 'time_of_day_24','scalar_axis'])
var5.name = 'areo'

var6 = xr.DataArray(DDS.taudust_VIS, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var6.name = 'taudust_VIS'

var7 = xr.DataArray(DS_fixed.zsurf, coords={'lon': DDS.lon, 'lat': DDS.lat},dims=['lat','lon'])
var7.name = 'zsurf'

NEW_DF = xr.merge([var1, var2, var3, var4,var5,var6,var7])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_savefile.nc')

# 2. Calculate Power [W or kW]  for Each Turbine Using Power Curve
# Enercon E33

print('Calculating E33 Power')
power333kw = (0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
thresh333kw = (0.,2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13.)
thresh333kw_top = (2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 28.)
tmp=test.veladj.sel(zagl=50)
pow_diurn_E33 = tmp
pow_diurn_E33 = np.where(tmp > 28, 0, pow_diurn_E33)              # first set cut off wind speed
for x,y,z in zip(thresh333kw,thresh333kw_top,power333kw):   # loop through
    print(x,y,z)
    pow_diurn_E33 = np.where(np.logical_and(tmp >=x, tmp< y),z,pow_diurn_E33)

# Aeolos-V
power300w = (0., 0.003, 0.005, 0.007, 0.014, 0.025, 0.043, 0.06, 0.08, 0.1, 0.12, 0.15, 0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.33, 0.36, 0.38, 0.4)
thresh300w = (0, 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5)
thresh300w_top = (1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 13.5)
tmp = test.veladj.sel(zagl=5)
pow_diurn_AAV = tmp
pow_diurn_AAV = np.where(tmp > 13, 0, pow_diurn_AAV)      # first set cut off wind speed
for x,y,z in zip(thresh300w, thresh300w_top, power300w):   # loop through 
    print(x,y,z)
    pow_diurn_AAV = np.where(np.logical_and(tmp >=x, tmp < y),z,pow_diurn_AAV)
        
# Jacobs 31-20
print('Calculating Jacobs 31-20 Power')
thresh20kw = (0.,3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,11.,11.5, 12., 12.5,13.,13.5)
thresh20kw_top = (3., 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,11.,11.5, 12., 12.5,13.,13.5,14.)
power20kw = (0.,0.09,0.26,0.5,0.83,1.22,1.88,2.6,3.43,4.5,5.6,6.79,8.0,9.65,11.45,13.1,15.,17.2,18.7,19.8,20,20,20.)
tmp=test.veladj.sel(zagl=30)
pow_diurn_Jac = tmp
pow_diurn_Jac = np.where(tmp > 14, 0, pow_diurn_Jac)              # first set cut off wind speed
for x,y,z in zip(thresh20kw,thresh20kw_top,power20kw):   # loop through
    print(x,y,z)
    pow_diurn_Jac = np.where(np.logical_and(tmp >=x, tmp < y),z,pow_diurn_Jac)

# NREL 5MW
print('Calculating NREL 5MW Power')
power5mw = (0.,40.5, 177.7, 403.9, 737.6, 1187.2, 1771.1, 2518.6, 3448.4, 4562.5, 5000.)
thresh5mw = (0.,3., 4., 5., 6., 7., 8., 9., 10., 11., 12.)
thresh5mw_top = (3.,4., 5., 6., 7., 8., 9., 10., 11., 25.)
tmp = test.veladj.sel(zagl=100)
pow_diurn_NREL = tmp
pow_diurn_NREL = np.where(tmp > 25, 0, pow_diurn_NREL)              # first set cut off wind speed
for x,y,z in zip(thresh5mw,thresh5mw_top,power5mw):   # loop through
    print(x,y,z)
    pow_diurn_NREL = np.where(np.logical_and(tmp >=x, tmp < y),z,pow_diurn_NREL)


# use atmospheric average file for short wave solar surface energy flux [W/m2]
var1 = xr.DataArray(DDS.areo, coords={'time': DDS.time,'time_of_day_24': DDS.time_of_day_24, 'scalar_axis': DDS.scalar_axis},dims=['time', 'time_of_day_24','scalar_axis'])
var1.name = 'areo'

var2 = xr.DataArray(pow_diurn_E33, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var2.name = 'pow_diurn_E33'

var3 = xr.DataArray(pow_diurn_Jac, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var3.name = 'pow_diurn_Jac'

var4 = xr.DataArray(pow_diurn_NREL, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var4.name = 'pow_diurn_NREL'

var5 = xr.DataArray(pow_diurn_AAV, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var5.name = 'pow_diurn_AAV'

NEW_DF = xr.merge([var1,var2,var3,var4,var5])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_pow2.nc')

# ========================================================
# Analysis of MY28
# ========================================================

dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/02094.fixed.nc'
DS_fixed = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY28_highres/*.atmos_diurn_T_zagl.min.nc'
DDS = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')

# 1. Calculate Wind Power Density
# calculate wind speed (again using diurnal file)
vel_diurn = sqrt(DDS.ucomp**2 + DDS.vcomp**2)

# calculate wind profile denisty [W/m2]
wpd = 0.5 * DDS.rho * vel_diurn**3


# use atmospheric average file for short wave solar surface energy flux [W/m2]
swflx = DDS.swflx

var1 = xr.DataArray(swflx, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var1.name = 'swflx'

var2 = xr.DataArray(wpd, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24, 'zagl': DDS.zagl},dims=['time','time_of_day_24','zagl','lat','lon'])
var2.name = 'wpd'

var3 = xr.DataArray(vel_diurn, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24, 'zagl': DDS.zagl},dims=['time', 'time_of_day_24','zagl','lat','lon'])
var3.name = 'ws'

var4 = xr.DataArray(DDS.areo, coords={'time': DDS.time,'time_of_day_24': DDS.time_of_day_24, 'scalar_axis': DDS.scalar_axis},dims=['time', 'time_of_day_24','scalar_axis'])
var4.name = 'areo'

var5 = xr.DataArray(DDS.taudust_VIS, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var5.name = 'taudust_VIS'

var6 = xr.DataArray(DS_fixed.zsurf, coords={'lon': DDS.lon, 'lat': DDS.lat},dims=['lat','lon'])
var6.name = 'zsurf'

NEW_DF = xr.merge([var1, var2, var3, var4,var5,var6])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY28_highres/MY28_highres_savefile.nc')

# 2. Calculate Power [W or kW]  for Each Turbine Using Power Curve
# Enercon E33
print('Calculating E33 Power')
power333kw = (0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
thresh333kw = (0.,2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13.)
thresh333kw_top = (2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 28.)
tmp=test.ws.sel(zagl=50)
pow_diurn_E33 = tmp
pow_diurn_E33 = np.where(tmp > 28, 0, pow_diurn_E33)              # first set cut off wind speed
for x,y,z in zip(thresh333kw,thresh333kw_top,power333kw):   # loop through
    print(x,y,z)
    pow_diurn_E33 = np.where(np.logical_and(tmp >=x, tmp< y),z,pow_diurn_E33)

# Jacobs 31-20
print('Calculating Jacobs 31-20 Power')
thresh20kw = (0.,3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,11.,11.5, 12., 12.5,13.,13.5)
thresh20kw_top = (3., 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5,11.,11.5, 12., 12.5,13.,13.5,14.)
power20kw = (0.,0.09,0.26,0.5,0.83,1.22,1.88,2.6,3.43,4.5,5.6,6.79,8.0,9.65,11.45,13.1,15.,17.2,18.7,19.8,20,20,20.)
tmp=test.ws.sel(zagl=30)
pow_diurn_Jac = tmp
pow_diurn_Jac = np.where(tmp > 14, 0, pow_diurn_Jac)              # first set cut off wind speed
for x,y,z in zip(thresh20kw,thresh20kw_top,power20kw):   # loop through
    print(x,y,z)
    pow_diurn_Jac = np.where(np.logical_and(tmp >=x, tmp < y),z,pow_diurn_Jac)

# NREL 5MW
print('Calculating NREL 5MW Power')
power5mw = (0.,40.5, 177.7, 403.9, 737.6, 1187.2, 1771.1, 2518.6, 3448.4, 4562.5, 5000.)
thresh5mw = (0.,3., 4., 5., 6., 7., 8., 9., 10., 11., 12.)
thresh5mw_top = (3.,4., 5., 6., 7., 8., 9., 10., 11., 25.)
tmp = test.ws.sel(zagl=100)
pow_diurn_NREL = tmp
pow_diurn_NREL = np.where(tmp > 25, 0, pow_diurn_NREL)              # first set cut off wind speed
for x,y,z in zip(thresh5mw,thresh5mw_top,power5mw):   # loop through
    print(x,y,z)
    pow_diurn_NREL = np.where(np.logical_and(tmp >=x, tmp < y),z,pow_diurn_NREL)


# use atmospheric average file for short wave solar surface energy flux [W/m2]
var1 = xr.DataArray(DDS.areo, coords={'time': DDS.time,'time_of_day_24': DDS.time_of_day_24, 'scalar_axis': DDS.scalar_axis},dims=['time', 'time_of_day_24','scalar_axis'])
var1.name = 'areo'

var2 = xr.DataArray(pow_diurn_E33, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var2.name = 'pow_diurn_E33'

var3 = xr.DataArray(pow_diurn_Jac, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var3.name = 'pow_diurn_Jac'

var4 = xr.DataArray(pow_diurn_NREL, coords={'lon': DDS.lon, 'lat': DDS.lat,'time': DDS.time, 'time_of_day_24': DDS.time_of_day_24},dims=['time', 'time_of_day_24','lat','lon'])
var4.name = 'pow_diurn_NREL'

NEW_DF = xr.merge([var1,var2,var3,var4])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY28_highres_pow.nc')

#===============================================================
# Calculate Percent Time Wind/Solar/Total Enercon E33 &
# 2500m2 Solar Array
# ==============================================================
# Import Base Datasets
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_pow2.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_savefile.nc'
DDS_swflx = xr.open_dataset(dataDIR, decode_times=False)

count_solar2500,count_total2500 = np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon))),np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon)))
for i in range(0,len(DDS_MY24.lat)):
    print(i)
    for j in range(0,len(DDS_MY24.lon)):
        count_solar2500[i,j] = np.count_nonzero(DDS_swflx.swflx[:,:,i,j]*0.2*2.5  >=24)
        count_total2500[i,j] = np.count_nonzero(DDS_swflx.swflx[:,:,i,j]*0.2*2.5 + DDS_MY24.pow_diurn_E33[:,:,i,j] >=24)


# Save variables
var1 = xr.DataArray(count_solar2500, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var1.name = 'count_solar2500'
var2 = xr.DataArray(count_total2500, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var2.name = 'count_total2500'

NEW_DF = xr.merge([var1,var2])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_solarwindpercent2500.nc')

#===============================================================
# Calculate Percent Time Wind/Solar/Total Alternative Turbine &
# 2500m2 Solar Array
# ==============================================================
#
count_NREL,count_Jac = np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon))),np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon)))
for i in range(0,len(DDS_MY24.lat)):
    print(i)
    for j in range(0,len(DDS_MY24.lon)):
        count_NREL[i,j] = np.count_nonzero(DDS_MY24.pow_diurn_NREL[:,:,i,j]+DDS_swflx.swflx[:,:,i,j]*0.2*2.5 >=24)
        count_Jac[i,j] = np.count_nonzero(DDS_MY24.pow_diurn_Jac[:,:,i,j]+DDS_swflx.swflx[:,:,i,j]*0.2*2.5>=24)

time = len(DDS_swflx.time)*len(DDS_swflx.time_of_day_24)
percent_NREL = count_NREL/time
percent_Jac = count_Jac/time

# Save variables
var1 = xr.DataArray(percent_Jac, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var1.name = 'percent_Jac'
var2 = xr.DataArray(percent_NREL, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var2.name = 'percent_NREL'

NEW_DF = xr.merge([var1,var2])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_percentaltturb_v2.nc')


#==============================================================
# Percent Time Wind/Solar Enercon E33 & 300m2 Solar Array
#==============================================================

count_solar300,count_total300 = np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon))),np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon)))
count_wind2 = np.zeros((len(DDS_MY24.lat),len(DDS_MY24.lon)))

for i in range(0,len(DDS_MY24.lat)):
    print(i)
    for j in range(0,len(DDS_MY24.lon)):
        count_wind2[i,j] = np.count_nonzero(DDS_MY24.pow_diurn_E33[:,:,i,j] >= 2.2)
        count_solar300[i,j] = np.count_nonzero(DDS_swflx.swflx[:,:,i,j]*0.2*0.3  >=24)
        count_total300[i,j] = np.count_nonzero(DDS_swflx.swflx[:,:,i,j]*0.2*0.3 + DDS_MY24.pow_diurn_E33[:,:,i,j] >=24)

time = len(DDS_swflx.time) * len(DDS_swflx.time_of_day_24)
percent_solar300 = count_solar300/time
percent_total300 = count_total300/time
percent_wind2 = count_wind2/time

# Save variables
var1 = xr.DataArray(percent_solar300, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var1.name = 'percent_solar300'
var2 = xr.DataArray(percent_total300, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var2.name = 'percent_total300'
var3 = xr.DataArray(percent_wind2, coords={'lon': DDS_MY24.lon, 'lat': DDS_MY24.lat},dims=['lat','lon'])
var3.name = 'percent_wind2'



NEW_DF = xr.merge([var1,var2,var3])
NEW_DF.to_netcdf('/Users/vhartwic/Documents/codes/VHplots/Data/fv3_2.0_windenergy/MY24_highres/MY24_highres_solarwindpercent300.nc')
