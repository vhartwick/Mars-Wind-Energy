# Figure S2 : Python script for base analysis and figures for Hartwick+2022
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
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_savefile.nc'
DDS_MY24 = xr.open_dataset(dataDIR, decode_times=False)

# -------------------------------------------------------------------------------
# FIGURE S2
# Annual average wind (top panels) and solar power density (bottom panels) [W/m2]
# in one-hour intervals of local time. All locations in each panel show the same
# local time simultaneously. Solar power exceeds 140 W/m2 during daytime hours
# (~7.5-15.5LT) while wind power maximizes at night (~21.5-8.5LT) and in the poles.
# -------------------------------------------------------------------------------
ct = 2
print('Making Figure S2')

# Find Annual Average
wpd_diurn_ave_MY24 = DDS_MY24.wpd.mean('time')
swflx_diurn_ave_MY24 = DDS_MY24.swflx.mean('time')

levels = np.arange(0,160,10)
levels_swflx = np.arange(0,500,25)

fig, axs = py.subplots(4,6, sharey=True, sharex=True, figsize=(10,8))
fig.subplots_adjust(wspace=0.05,hspace=0.35)
fig.suptitle('(a)', ha='left')

for i in range(6):
        cs = axs[0,i].contourf(DDS_MY24.lon,DDS_MY24.lat,wpd_diurn_ave_MY24.sel(zagl=50,time_of_day_24=i+0.5),levels=levels,cmap=py.cm.viridis)
        axs[0,i].set_title(str(i+0.5)+'LT',fontsize=10)
        if i == 0:
                axs[0,i].set_ylabel('Latitude',fontsize=8)
                axs[0,i].xaxis.set_major_locator(MultipleLocator(60))
                axs[0,i].yaxis.set_major_locator(MultipleLocator(30))
                axs[0,i].set_yticks(-90,-60,-30,0,30,60,90)
                axs[0,i].set_yticklabels(['','-60','','0','','60',''])
for i in range(6):
        cs = axs[1,i].contourf(DDS_MY24.lon,DDS_MY24.lat,wpd_diurn_ave_MY24.sel(zagl=50,time_of_day_24=i+6.5),levels=levels,cmap=py.cm.viridis)
        axs[1,i].set_title(str(i+6.5)+'LT',fontsize=10)
        if i == 0:
                axs[1,i].set_ylabel('Latitude',fontsize=8)
for i in range(6):
        cs = axs[2,i].contourf(DDS_MY24.lon,DDS_MY24.lat,wpd_diurn_ave_MY24.sel(zagl=50,time_of_day_24=i+12.5),levels=levels,cmap=py.cm.viridis)
        axs[2,i].set_title(str(i+12.5)+'LT',fontsize=10)
        if i == 0:
                axs[2,i].set_ylabel('Latitude',fontsize=8)

for i in range(6):
        cs = axs[3,i].contourf(DDS_MY24.lon,DDS_MY24.lat,wpd_diurn_ave_MY24.sel(zagl=50,time_of_day_24=i+18.5),levels=levels,cmap=py.cm.viridis)
        axs[3,i].set_title(str(i+18.5)+'LT',fontsize=10)
        axs[3,i].set_xlabel('Longitude',fontsize=8)
        axs[3,i].set_xticks(0,60,120,180,240,300,360)
        axs[3,i].set_xticklabels(['','60','','180','','300',''])
        if i == 0:
                axs[3,i].set_ylabel('Latitude',fontsize=8)


norm= matplotlib.colors.Normalize(vmin=cs.cvalues.min(), vmax=cs.cvalues.max())
sm = py.cm.ScalarMappable(norm=norm, cmap = cs.cmap)
sm.set_array([])

# add axis, (left, bottom, width, height)
cbar_ax = py.gcf().add_axes([0.93, 0.11, 0.01, 0.77])
clb = fig.colorbar(sm,  cax=cbar_ax, orientation='vertical')
clb.set_label('$\ [W/m^{2}]$')
py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

# Now Solar
print('Now, Solar')

fig, axs = py.subplots(4,6, sharey=True, sharex=True, figsize=(10,8))
fig.subplots_adjust(wspace=0.05,hspace=0.35)
fig.suptitle('(b)', ha='left')

for i in range(6):
        cs = axs[0,i].contourf(DDS_MY24.lon,DDS_MY24.lat,swflx_diurn_ave_MY24.sel(time_of_day_24=i+0.5),levels=levels_swflx,cmap=py.cm.viridis)
        axs[0,i].set_title(str(i+0.5)+'LT',fontsize=10)
        if i == 0:
                axs[0,i].set_ylabel('Latitude',fontsize=8)
                axs[0,i].xaxis.set_major_locator(MultipleLocator(60))
                axs[0,i].yaxis.set_major_locator(MultipleLocator(30))
                axs[0,i].set_yticklabels(['','-60','','0','','60',''])
for i in range(6):
        cs = axs[1,i].contourf(DDS_MY24.lon,DDS_MY24.lat,swflx_diurn_ave_MY24.sel(time_of_day_24=i+6.5),levels=levels_swflx,cmap=py.cm.viridis)
        axs[1,i].set_title(str(i+6.5)+'LT',fontsize=10)
        if i == 0:
                axs[1,i].set_ylabel('Latitude',fontsize=8)
for i in range(6):
        cs = axs[2,i].contourf(DDS_MY24.lon,DDS_MY24.lat,swflx_diurn_ave_MY24.sel(time_of_day_24=i+12.5),levels=levels_swflx,cmap=py.cm.viridis)
        axs[2,i].set_title(str(i+12.5)+'LT',fontsize=10)
        if i == 0:
                axs[2,i].set_ylabel('Latitude',fontsize=8)

for i in range(6):
        cs = axs[3,i].contourf(DDS_MY24.lon,DDS_MY24.lat,swflx_diurn_ave_MY24.sel(time_of_day_24=i+18.5),levels=levels_swflx,cmap=py.cm.viridis)
        axs[3,i].set_title(str(i+18.5)+'LT',fontsize=10)
        axs[3,i].set_xlabel('Longitude',fontsize=8)
        axs[3,i].set_xticklabels(['','60','','180','','300',''])
        if i == 0:
                axs[3,i].set_ylabel('Latitude',fontsize=8)


norm= matplotlib.colors.Normalize(vmin=cs.cvalues.min(), vmax=cs.cvalues.max())
sm = py.cm.ScalarMappable(norm=norm, cmap = cs.cmap)
sm.set_array([])
# add axis, (left, bottom, width, height)
cbar_ax = py.gcf().add_axes([0.93, 0.11, 0.01, 0.77])

clb = fig.colorbar(sm,  cax=cbar_ax, orientation='vertical')
clb.set_label('$\ [W/m^{2}]$')
py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}b.eps',dpi=300)

# Send out Data
var1 = xr.DataArray(swflx_diurn_ave_MY24, coords={'lon':DDS_MY24.lon,'lat':DDS_MY24.lat,'time_of_day_24':DDS_MY24.time_of_day_24},dims=['time_of_day_24','lat','lon'])
var1.name = 'annual_ave_pow_E33'

var2 = xr.DataArray(wpd_diurn_ave_MY24, coords={'lon':DDS_MY24.lon,'lat':DDS_MY24.lat,'zagl':DDS_MY24.zagl,'time_of_day_24':DDS_MY24.time_of_day_24},dims=['time_of_day_24','zagl','lat','lon'])
var2.name = 'wpd_diurn_ave_MY24'

NEW_DF = xr.merge([var1,var2])
NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')
