# Figure S4: Python script for base analysis and figures for Hartwick+2022
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
DDS_swflx = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S4
# Wind turbines generate significant energy on an annual basis. Solar panel
# arrays would need to exceed 600 m2 on average to match energy produced by one
# turbine. (a) Annual energy production (GWh) by the Enercon E33 wind turbine
# over one Mars year compared with (b) the solar panel array dimensions (m2)
# needed to generate an equivalent amount of energy.
# -------------------------------------------------------------------------------
ct = 4
print('Making Figure S4')

# Calculate Annual Energy Production
pow_tot_diurn_MY24 = np.sum(DDS_MY24.pow_diurn_E33, axis=1) # sum over time of day
energy_per_sol = pow_tot_diurn_MY24.mean('time')
pow_tot_diurn_MY24 = np.sum(pow_tot_diurn_MY24, axis=0)
aep_MY24 = pow_tot_diurn_MY24*5.


# Print AEP at Highlighted ROI for Table S3
lat_roi = [35,29,38,-23,-31,-11,-36,-3,-16,-28]
lon_roi = [23,17,48,78,94,174,200,325,311,97]

print('max',aep_MY24.max(),aep_MY24.argmax(...))
print('max', energy_per_sol.max(),energy_per_sol.argmax(...))
print('lat=', DDS_MY24.lat[115], 'lon=',DDS_MY24.lon[137])
weights = np.cos(np.deg2rad(DDS_MY24.lat))
weights.name = "weights"
print('global average aep=', aep_MY24.weighted(weights).mean(('lon','lat')))
#for i in range(0,len(lat_roi)):
#	print(aep_MY24.sel(lat=lat_roi[i],lon=lon_roi[i],method='nearest'))
#	print('per_sol=',energy_per_sol.sel(lat=lat_roi[i],lon=lon_roi[i],method='nearest'))

# Calculate Equivalent Solar Panel Array Area
pow_m2_solar = DDS_swflx.swflx*0.2 * (1/1000)   # W to kW
spow_tot_diurn_MY24 = np.sum(pow_m2_solar, axis=1)
senergy_per_sol = spow_tot_diurn_MY24.mean('time')
spow_tot_diurn_MY24 = np.sum(spow_tot_diurn_MY24, axis=0)
saep_MY24 = spow_tot_diurn_MY24 * 5 # no 1000 bc swflx units are in W not kW
S_solar = aep_MY24/saep_MY24
print('S_solar max=', S_solar.max(),S_solar.argmax(...))
print('senergy max=', spow_tot_diurn_MY24.max(),spow_tot_diurn_MY24.argmax(...))
print('senergy_per_sol max=', senergy_per_sol.max(),senergy_per_sol.argmax(...))
print('lat=',DDS_MY24.lat[92],'lon=',DDS_MY24.lon[67])
print('global average saep=', saep_MY24.weighted(weights).mean(('lon','lat')))
# Plot
energylim=[0,15,24,35,100]  # energy requirement per sol (Rucker 2015 https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20150000526.pdf)
slevels = np.arange(0,1600,100)
#slevels = np.arange(0,1000,50)
smin = np.arange(0,1600,400)
fig, ax = py.subplots(1,2, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)

fig1 = ax[0].contourf(DDS_MY24.lon, DDS_MY24.lat, aep_MY24[:,:]*1e-6, 20, cmap=py.cm.viridis)
ax[0].set_title('(a) AEP [GWh]')
ax[0].set_ylabel('Latitude')
ax[0].set_xlabel('Longitude')
ax[0].xaxis.set_major_locator(MultipleLocator(60))
ax[0].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig1,ax=ax[0])
cbar.set_label('[GWh]')

cmap = matplotlib.cm.plasma
norm = matplotlib.colors.BoundaryNorm(slevels,cmap.N,extend='max')
fig2 = ax[1].contourf(DDS_MY24.lon, DDS_MY24.lat, S_solar[:,:], levels=slevels,cmap=py.cm.magma, norm=norm, extend='max')
ax[1].set_title('(b) Solar Array Size $\ [m^{2}$]')
cbar = py.colorbar(fig2,extend='max',ax=ax[1])
cbar.set_label('$\ [m^{2}$]')
fig2 = ax[1].contour(DDS_MY24.lon, DDS_MY24.lat, S_solar[:,:], levels=smin, cmap=py.cm.Greys)
py.clabel(fig2, fontsize=9,inline=1, fmt='%1.0f')
ax[1].set_xlabel('Longitude')

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

#Send out Data
var1 = xr.DataArray(aep_MY24, coords={'lon':DDS_MY24.lon, 'lat':DDS_MY24.lat},dims=['lat','lon'])
var1.name = 'aep_MY24'

var2 = xr.DataArray(S_solar, coords={'lon':DDS_MY24.lon,'lat':DDS_MY24.lat},dims=['lat','lon'])
var2.name = 'S_solar'

NEW_DF = xr.merge([var1,var2])
NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')

