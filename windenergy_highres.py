# Create wind energy plots from HIGH RES simulations
# Victoria Hartwick
# Updated 12/17/21

# Import Modules
import os
import matplotlib
import matplotlib.pyplot as py
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.cm import get_cmap   #contour color tables
from matplotlib.ticker import MultipleLocator, FuncFormatter #format ticks
from netCDF4 import Dataset, MFDataset
from numpy import sqrt, exp, max, mean, min, log, log10,int
#from bokeh.io import output_notebook
#output_notebook()  # set to output the plot in the notebook
#from bokeh.plotting import show, figure
from itertools import product
import xarray as xr
from matplotlib.patches import Rectangle
from dbfread import DBF

# Get Current Path
PATH = os.getcwd()
print(f'\nCurrent path is:\n{PATH}\n')

# Open first set of files
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/*.atmos_diurn_T_zagl_min.nc'
DDS_SY1 = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')
   
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/02094.fixed.nc'
DS_fixed = xr.open_dataset(dataDIR, decode_times=False)

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
    idx = (np.abs(DDS_SY1.lon-tmp[i])).argmin()
    lon_idx += [idx]
    
tmp = np.array(latitude)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_SY1.lat-tmp[i])).argmin()
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


# calculate wind speed (again using diurnal file)
vel_diurn_SY1 = sqrt(DDS_SY1.ucomp**2 + DDS_SY1.vcomp**2)

# calculate wind profile denisty [W/m2]
wpd_SY1 = 0.5 * DDS_SY1.rho * vel_diurn_SY1**3

# use atmospheric average file for short wave solar surface energy flux [W/m2]
swflx_diurn_SY1 = DDS_SY1.swflx

# calculate annual and diurnal averages
tmp = wpd_SY1.mean('time')
wpd_annual_ave_SY1 = tmp.mean('time_of_day_24')
tmp = swflx_diurn_SY1.mean('time')
swflx_annual_ave_SY1 = tmp.mean('time_of_day_24')
taudust_annual_ave_SY1 = DDS_SY1.taudust_VIS.mean('time')

wpd_time_ave_SY1 = wpd_SY1.mean('time_of_day_24')
wpd_diurn_ave_SY1 = wpd_SY1.mean('time')
swflx_diurn_ave_SY1 = swflx_diurn_SY1.mean('time')

# --------------------------------------------------------------------------------
# FIGURE 1
# --------------------------------------------------------------------------------
print('Making Figure 1')
ct = 1
levels = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70]
levels_swf = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]

fig, ([ax1,ax2,ax3]) = py.subplots(1,3, sharey=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)

fig1 = ax1.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_annual_ave_SY1.sel(zagl=5), cmap=py.cm.viridis, levels=levels_swf)
ax1.contour(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf,5, cmap=py.cm.Greys_r)
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
clb=py.colorbar(fig1, ax=ax1, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

fig2 = ax2.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_annual_ave_SY1.sel(zagl=30), cmap=py.cm.viridis, levels=levels_swf)
ax2.contour(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf,5, cmap=py.cm.Greys_r)
ax2.set_xlabel('Longitude')
clb=py.colorbar(fig2, ax=ax2, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

fig3 = ax3.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_annual_ave_SY1.sel(zagl=100), cmap=py.cm.viridis, levels=levels_swf)
ax3.contour(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf,5, cmap=py.cm.Greys_r)
ax3.set_xlabel('Longitude')
clb=py.colorbar(fig3, ax=ax3, orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

py.savefig(f'{PATH}/windenergy_fig{ct}a.eps',dpi=200)
# 50m
fig, ([ax1,ax2]) = py.subplots(1,2, sharey=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)
fig1 = ax1.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_annual_ave_SY1.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax1.contour(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf,5, cmap=py.cm.Greys_r)
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
clb = py.colorbar(fig1,ax=ax1)
clb.set_label('$\ [W/m^{2}$]')

# compare with Solar Energy
fig2 = ax2.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_annual_ave_SY1.sel(zagl=50)/swflx_annual_ave_SY1,cmap=py.cm.viridis)
ax2.contour(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf,5, cmap=py.cm.Greys_r)
ax2.set_xlabel('Longitude')
clb = py.colorbar(fig2, ax=ax2)

py.savefig(f'{PATH}/windenergy_fig{ct}b.eps',dpi=200)

# --------------------------------------------------------------------------------
# FIGURE 2
# --------------------------------------------------------------------------------
print('Making Figure 2')
ct = 2
# note number of years
yrs = [0,1,2,3]

# locate cardinal seasons in Ls array
areo = DDS_SY1.areo[:,0] 
#areo_SY1 = Dareo - 360*3
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

# average all years to produce multiyear seasonal averages (MC_season)
wpd_season_MY24 = []
swflx_season_MY24 = []
taudust_season_MY24 = []

tmp = wpd_SY1
tmp2 = np.mean(tmp[ls90,:,:,:],axis=0)
tmp3 = np.mean(tmp[ls180,:,:,:],axis=0)
tmp4 = np.mean(tmp[ls270,:,:,:],axis=0)
tmp5 = np.mean(tmp[ls0,:,:,:],axis=0)

tmp2 = tmp2.mean('time_of_day_24')
tmp3 = tmp3.mean('time_of_day_24')
tmp4 = tmp4.mean('time_of_day_24')
tmp5 = tmp5.mean('time_of_day_24')

wpd_season_MY24 += [tmp2.rename('wpd_Ls90')]
wpd_season_MY24 += [tmp3.rename('wpd_Ls180')]
wpd_season_MY24 += [tmp4.rename('wpd_Ls270')]
wpd_season_MY24 += [tmp5.rename('wpd_Ls0')]
wpd_MY24 = xr.merge(wpd_season_MY24)

# SOLAR
tmp = DDS_SY1.swflx
tmp2 = np.mean(tmp[ls90,:,:],axis=0)
tmp3 = np.mean(tmp[ls180,:,:],axis=0)
tmp4 = np.mean(tmp[ls270,:,:],axis=0)
tmp5 = np.mean(tmp[ls0,:,:],axis=0)

tmp2 = tmp2.mean('time_of_day_24')
tmp3 = tmp3.mean('time_of_day_24')
tmp4 = tmp4.mean('time_of_day_24')
tmp5 = tmp5.mean('time_of_day_24')

swflx_season_MY24 += [tmp2.rename('swflx_Ls90')]
swflx_season_MY24 += [tmp3.rename('swflx_Ls180')]
swflx_season_MY24 += [tmp4.rename('swflx_Ls270')]
swflx_season_MY24 += [tmp5.rename('swflx_Ls0')]
swflx_MY24 = xr.merge(swflx_season_MY24)

# DUST OPTICAL DEPTH
tmp = DDS_SY1.taudust_VIS
tmp2 = np.mean(tmp[ls90,:,:],axis=0)
tmp3 = np.mean(tmp[ls180,:,:],axis=0)
tmp4 = np.mean(tmp[ls270,:,:],axis=0)
tmp5 = np.mean(tmp[ls0,:,:],axis=0)
taudust_season_MY24 += [tmp2.rename('taudust_Ls90')]
taudust_season_MY24 += [tmp3.rename('taudust_Ls180')]
taudust_season_MY24 += [tmp4.rename('taudust_Ls270')]
taudust_season_MY24 += [tmp5.rename('taudust_Ls0')]
taudust_MY24 = xr.merge(taudust_season_MY24)

# plot
# Cardinal Seasons FV3 (Using DIURN FILE)
import warnings; warnings.simplefilter('ignore')

levels_swf = np.arange(0,160,10)
levels_swf2 = np.arange(0,240,20)
level_lim = [100,500]
dust_lim = [0,0.5] #is there a way to calculate this? optical depth when solar cell starts working at half efficiency or something?

# Plot WPD
fig, ([ax1,ax2],[ax3,ax4]) = py.subplots(2,2, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)
#fig.suptitle('Seasonal Average Wind Power Density at 50m $\ [W/m^{2}]$, MY 24')

fig1 = ax1.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_MY24.wpd_Ls90.sel(zagl=50).squeeze(), levels=levels_swf,cmap=py.cm.viridis)
ax1.set_title('Ls=80-100')
ax1.set_ylabel('Latitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))
clb = py.colorbar(fig1, ax=ax1)
clb.set_label('$\ [W/m^{2}$]')

fig1 = ax1.contour(DDS_SY1.lon, DDS_SY1.lat, swflx_MY24.swflx_Ls90, levels=levels_swf2, cmap=py.cm.plasma)
py.clabel(fig1, fontsize=9,inline=1, fontcolor='black', fmt='%1.0f')


fig2 = ax2.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_MY24.wpd_Ls180.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax2.set_title('Ls=170-190')
ax2.xaxis.set_major_locator(MultipleLocator(60))
ax2.yaxis.set_major_locator(MultipleLocator(30))
clb= py.colorbar(fig2, ax=ax2)
clb.set_label('$\ [W/m^{2}$]')

fig2 = ax2.contour(DDS_SY1.lon, DDS_SY1.lat, swflx_MY24.swflx_Ls180, levels=levels_swf2, cmap=py.cm.plasma)
py.clabel(fig2, fontsize=9,inline=1, fontcolor='black', fmt='%1.0f')


fig3 = ax3.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_MY24.wpd_Ls270.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax3.set_title('Ls=260-280')
ax3.set_ylabel('Latitude')
ax3.set_xlabel('Longitude')
ax3.xaxis.set_major_locator(MultipleLocator(60))
ax3.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig3, ax=ax3)
clb.set_label('$\ [W/m^{2}$]')

fig3 = ax3.contour(DDS_SY1.lon, DDS_SY1.lat, swflx_MY24.swflx_Ls270, levels=levels_swf2, cmap=py.cm.plasma)
py.clabel(fig3, fontsize=9,inline=1, fontcolor='black', fmt='%1.0f')


fig4 = ax4.contourf(DDS_SY1.lon, DDS_SY1.lat, wpd_MY24.wpd_Ls0.sel(zagl=50),levels=levels_swf,cmap=py.cm.viridis)
ax4.set_title('Ls=350-0')
ax4.set_xlabel('Longitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
clb=py.colorbar(fig4, ax=ax4)
clb.set_label('$\ [W/m^{2}$]')

fig4 = ax4.contour(DDS_SY1.lon, DDS_SY1.lat, swflx_MY24.swflx_Ls0, levels=levels_swf2, cmap=py.cm.plasma)
py.clabel(fig4, fontsize=9,inline=1,fontcolor='black', fmt='%1.0f')

py.savefig(f'{PATH}/windenergy_fig{ct}.eps',dpi=200)

#-------------------------------------------------------------------------------
# FIGURE 4
#-------------------------------------------------------------------------------

print('Making Figure 3')
ct = 4

# genesate day and night masks
m = np.ma.masked_where(DDS_SY1.swflx>0,wpd_SY1.sel(zagl=50))
print(m.shape)
night = np.mean(m,axis=0)
night = np.mean(night,axis=0)
print(night.shape)

m = np.ma.masked_where(DDS_SY1.swflx==0,wpd_SY1.sel(zagl=50))
day = np.mean(m,axis=0)
day = np.mean(day,axis=0)

# find nighttime ave
levels=np.arange(0,220,20)
levels_diff = [-300,-250,-200,-150,-100,-50,0,50,100,150,200,250,300]
levels_diff = [-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
levels_diff = [-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180,200]

fig, axs = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.1)

#fig.suptitle('Annual Average Nighttime (17-6 LT) Power [W/m2]')
atm_lev=50
axs = axs.ravel()


fig1 = axs[0].contourf(DDS_SY1.lon,DDS_SY1.lat,night,levels=levels,cmap=py.cm.viridis)
axs[0].set_title('Nighttime Wind')
axs[0].set_xlabel('Longitude')
axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[0].yaxis.set_major_locator(MultipleLocator(30))
   
#add colorabar
clb=py.colorbar(fig1,ax=axs[0],orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')
axs[0].set_ylabel('Latitude')

fig = axs[1].contourf(DDS_SY1.lon,DDS_SY1.lat,day,levels=levels,cmap=py.cm.viridis)
axs[1].set_title('Daytime Wind')
axs[1].set_xlabel('Longitude')
clb=py.colorbar(fig,ax=axs[1],orientation='horizontal')
clb.set_label('$\ [W/m^{2}$]')

# diff
fig = axs[2].contourf(DDS_SY1.lon,DDS_SY1.lat,night/day,cmap=py.cm.plasma,levels=np.arange(0,5.5,0.5))
axs[2].set_title('Night/Day')
axs[2].set_xlabel('Longitude')
py.colorbar(fig,ax=axs[2],orientation='horizontal')

py.savefig(f'{PATH}/windenergy_fig{ct}.eps',dpi=200)

#----------------------------------------------------------------------------------
# FIGURE 5
#----------------------------------------------------------------------------------

print('Making Figure 5')
ct = 5
dataDIR = '/Users/vhartwic/Downloads/Enercon_330kW_power_MY24_HIGHRES.nc'
DDS_E33 = xr.open_dataset(dataDIR, decode_times=False)

pow_m2_solar = DDS_SY1.swflx*0.2
spow_tot_diurn_MY24 = np.sum(pow_m2_solar, axis=0)
spow_tot_diurn_MY24 = np.sum(spow_tot_diurn_MY24, axis=0)
saep_MY24 = spow_tot_diurn_MY24*3600*5*3.6e-6 # no 1000 bc swflx units are in W not kW

 # AEP -  (from total)
pow_tot_diurn_MY24 = np.sum(DDS_E33.pow_diurn_SY1, axis=0)
pow_tot_diurn_MY24 = np.sum(pow_tot_diurn_MY24, axis=0)
aep_MY24 = pow_tot_diurn_MY24*1000*3600*5*3.6e-6 # try assumuming it is really units of kW

# annual average energy [kW] & AEP calculated from AAE for ref/diagnostics (want to be close to or identical to calc form total)
pow_ave_MY24 = np.mean(DDS_E33.pow_diurn_SY1, axis=0)
aae_MY24 = np.mean(pow_ave_MY24, axis=0)
aep_ave_MY24 = aae_MY24 * 1000 * 88775 * 5 * len(DDS_E33.time) * 3.6e-6 # average x sec/sol x 5sol/step * #steps * 1000 W/kw * 3.6e-6 kWh/J

# calculate array size to genesrate same AEP as wind
# solar aep = saep*area
#if we want solar aep to equal wind aep then awp_MY24 = saep_MY23 * S
#S = aep_MY24/saep_MY24
S_solar = aep_MY24/saep_MY24

# also calculate annual average power for E333 and general turbine, plus the energy/sol
energylim=[0,15,24,35,100]  # energy requirement per sol (Rucker 2015 https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20150000526.pdf)
slevels = np.arange(0,1600,100)
smin = np.arange(0,1600,400)
#smin = [1,10,100,1000]
fig, ax = py.subplots(1,2, sharey=True, sharex=True, figsize=(12,5))
#fig, ax = py.subplots(2,1, sharey=True, sharex=True, figsize=(8,10))
fig.subplots_adjust(wspace=0.1)

    
# PLOT
fig1 = ax[0].contourf(DDS_E33.lon, DDS_E33.lat, aep_MY24[:,:]*1e-6, 20, cmap=py.cm.viridis)
ax[0].set_title('AEP [GWh]',fontsize=16)
ax[0].set_ylabel('Latitude',fontsize=14)
ax[0].set_xlabel('Longitude',fontsize=14)

ax[0].xaxis.set_major_locator(MultipleLocator(60))
ax[0].yaxis.set_major_locator(MultipleLocator(30))
cbar=py.colorbar(fig1,ax=ax[0])
cbar.set_label('[GWh]',fontsize=14)


fig2 = ax[1].contourf(DDS_E33.lon, DDS_E33.lat, S_solar[:,:], levels=slevels, cmap=py.cm.magma)
ax[1].set_title('Solar Array Size $\ [m^{2}$]',fontsize=16)
cbar = py.colorbar(fig2,ax= ax[1])
cbar.set_label('$\ [m^{2}$]',fontsize=14)
fig2 = ax[1].contour(DDS_E33.lon, DDS_E33.lat, S_solar[:,:], levels=smin, cmap=py.cm.Greys)
py.clabel(fig2, fontsize=9,inline=1,fontcolor='black', fmt='%1.0f')
#ax[1].set_ylabel('Latitude',fontsize=14)
ax[1].set_xlabel('Longitude',fontsize=14)

# First we need to convert the latitudes to radians
latr = np.deg2rad(DDS_SY1.lat)

# Use the cosine of the converted latitudes as weights for the average
weights = np.cos(latr)

print(np.max(S_solar).values)
tmp = S_solar.mean('lon')
tmp2 = np.average(tmp,weights=weights)
print(tmp2)

py.savefig(f'{PATH}/windenergy_fig{ct}.eps',dpi=200)

#-------------------------------------------------------------------------------
# FIGURE 6
#-------------------------------------------------------------------------------

print('Making Figure 6')
ct = 6

# plot all locations with diurnal average power greater than some limit (ALL THE TIME)
# choose locations with energy (kW) always greater than the minimum energy requirements (24) 
# or the average greater than the mim energy requirement

# diurnal mean
tmp_SY1 = np.mean(DDS_E33.pow_diurn_SY1,axis=1) 


lat_sav24, lon_sav24 = [], []
lat_sav15, lon_sav15 = [], []
lat_sav2, lon_sav2 = [], []

for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp = tmp_SY1[:,i,j]
        if np.any(tmp < 24):
            continue
        else:
            lat_sav24+= [DDS_E33.lat[i]]
            lon_sav24 +=[DDS_E33.lon[j]]

for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp = tmp_SY1[:,i,j]
        if np.any(tmp < 14.2):
            continue
        else:
            lat_sav15 += [DDS_E33.lat[i]]
            lon_sav15 +=[DDS_E33.lon[j]]

for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp = tmp_SY1[:,i,j]
        if np.any(tmp < 2.2):
            continue
        else:
            lat_sav2 += [DDS_E33.lat[i]]
            lon_sav2 +=[DDS_E33.lon[j]]

# Regions that always produce daily average power greater than some value
# highly restrictive so also plot annual average power > 24
tmp = DDS_E33.pow_diurn_SY1.mean('time')
annual_ave_pow_E33 = tmp.mean('time_of_day_24')
level_lim = [0,14.2,24,1000]

fig, ax4 = py.subplots(1,1,sharey=True, figsize=(15,8))
fig4 = ax4.contourf(DDS_E33.lon, DDS_E33.lat, DS_fixed.zsurf, 20,cmap=py.cm.viridis)
ax4.plot(lon_sav2, lat_sav2, 'ok', markersize=10,fillstyle='none', label='>2.2kW')
ax4.plot(lon_sav15, lat_sav15, 'vk', markersize=15, markerfacecolor='blue',label='>14.2kW')
ax4.plot(lon_sav24, lat_sav24, 'sk', markersize=15, markerfacecolor='red', label='>24kW')
ax4.contour(DDS_E33.lon, DDS_E33.lat, annual_ave_pow_E33, cmap=py.cm.plasma_r, levels=level_lim)
#ax4.set_title('Proposed Landing Sites')
ax4.set_xlabel('Longitude')
ax4.set_ylabel('Latitude')
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))


# add landing sites
ax4.plot(longitude, latitude, 'ok',label='ROIs')
#ax4.plot(nhlim_lon,nhlim_lat,'--k')
#ax4.plot(shlim_lon,shlim_lat,'--k')
ax4.legend(loc='upper left')

# labels
#ax4.text(247, 47, 'Melon+2004',color='black', fontsize=13)
#ax4.text(20, -55, 'Melon+2004',color='black', fontsize=13)

py.savefig(f'{PATH}/windenergy_fig{ct}.eps,dpi=200')

#----------------------------------------------------------------------------
# FIGURE 7
#----------------------------------------------------------------------------

print('Making Figure 7')
ct=7

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


#Melon+2004 (NH)
nhlim_lon = [2.607548276295404, 15.019478071461492,38.31148159008086, 53.86629246429021, 67.40166715536378,64.74636618774349, 
73.9651489129979, 84.97591421270891, 98.20215305994219, 107.89343609935912, 110.87420935785195, 
119.8165291333305, 124.08411175805298, 137.37569639341518, 160.42265320655133, 160.60235412390568, 168.51799103589826, 
174.00703723872152,180.42253521126761, 181.6901408450704, 189.29577464788733,199.43661971830986, 210, 217.6056338028169, 224.7887323943662, 
237.88732394366195, 240,248.45070422535213, 261.9718309859155, 272.11267605633805, 275.49295774647885, 281.40845070422534, 
283.5211267605634,292.8169014084507, 302.11267605633805, 303.38028169014086, 318.59154929577466, 330, 338.87323943661977,
351.9718309859155, 358]
nhlim_lat=[47.189209567293624,44.506304193021435,44.55007749340262, 46.94726259791399, 47.3109789301722, 49.67394127256737
, 51.72093997402927, 51.741632806936735, 49.73681564947852, 50.76986553847446, 50.43718845557746, 49.43915720688647
, 48.770619528337484, 45.41280945000632, 50.53030620366108, 46.809575671260454, 49.530683198592584, 49.87927784526456
,  48.568047337278095, 48.23076923076922, 46.20710059171596, 44.520710059171584, 44.183431952662716, 42.834319526627205
, 38.449704142011825, 44.520710059171584, 44.85798816568046, 45.53254437869821, 42.15976331360945, 42.15976331360945
, 46.881656804733716, 46.881656804733716, 50.25443786982247, 47.21893491124259, 48.90532544378697, 52.6153846153846
, 52.6153846153846, 52.27810650887572, 51.94082840236685, 49.57988165680472, 48.568047337278095]

#Melon+2004 (SH)
shlim_lon=[2, 4.624471159887719, 14.039291249528759, 30.539102752062973,37.673103506052854,48.58585012357054, 57.34846898169479, 
72.18698948603026, 73.89728982532569, 84.4845641519708, 92.62891132241441, 107.35307669752439, 130.22158924307791, 
137.1117999413563, 146.41226490177183, 157.4720395425794, 169.36245968248647, 179.1190885100322,
183.38028169014083, 209.57746478873239, 221.83098591549296, 226.90140845070422, 236.19718309859155, 240.42253521126761, 
247.6056338028169, 261.1267605633803, 272.53521126760563, 282.2535211267606, 293.2394366197183, 319.85915492957747, 
329.1549295774648, 340.9859154929577, 351.1267605633803, 358]
shlim_lat=[-47.19637247099233, -47.18682193272734, -49.19880199388429, -48.82951451430485, -47.46299166422317, -45.41280945000628
, -42.690110166296634, -43.000502659908676, -43.3355673773719, -43.315670422653184, -45.33003811837641, -43.27269300046075
, -43.22971557826831, -45.584719138775995, -45.22896158840531, -46.223013446152535, -46.87722531730405, -47.19716834918104
, -48.5680473372781, -48.90532544378701, -48.5680473372781, -48.5680473372781, -48.230769230769226, -48.5680473372781
, -49.242603550295854, -49.242603550295854, -47.89349112426035, -47.89349112426035, -48.5680473372781, -48.5680473372781
, -46.544378698224875, -46.20710059171597, -47.556213017751475, -47.556213017751475]

# crater locations
lon_crater = [55.4631379962193, 76.89981096408323, 137.80718336483935, 175.9168241965974,
            194.07272727272726,199.96363636363637,190.8, 323.0181818181818, 
            337.41818181818184, 351.81818181818187] 

lat_crater = [ -11.121951219512184, 14.18048780487807,-3.6146341463414444, -12.234146341463415,
            -24.962529274004694, -36.17564402810305, -46.32084309133489, -31.370023419203747, 
            18.555035128805613, -0.13348946135830886]

# define ROIs
lon_new = [22.117202268431015, 22.797731568998103, 17.0132325141777, 47.97731568998114, 
             77.58034026465029, 96.63516068052931, 94.25330812854443, 
             102.07939508506618,104.12098298676747, 137.80718336483935,
             158.22306238185254, 176.93761814744803,174.21550094517963,175.9168241965974,
             199.96363636363637,220.25454545454545,300.1090909090909,
             313.8545454545455, 325.30909090909097,310.9090909090909] 

lat_new = [36.702439024390266, 34.75609756097563, 29.195121951219527, 37.81463414634148, 
            -23.07804878048779, -27.526829268292673, -31.41951219512194,
            -34.756097560975604, -35.868292682926835, -3.6146341463414444,
            -0.834146341463395, -10.009756097560967, -10.843902439024376, -12.234146341463415,
            -36.17564402810305, 35.64168618266979, -25.229508196721312, 9.744730679156902,
            -3.3372365339578494,-16.4192037470726]

# set contours & color
energylim = [0,14,24,100]
colors = ['white', 'grey', py.cm.viridis(0.75), py.cm.viridis(1)]

fig, (ax4) = py.subplots(nrows=1, figsize=(15,8))
ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
datafile = '/Users/vhartwic/Desktop/SWIM2.png'
fig=ax4.contourf(DDS_E33.lon, DDS_E33.lat, aae_MY24[:,:], cmap=py.cm.Greys, levels=energylim, alpha=.5)
fig=ax4.contour(DDS_E33.lon, DDS_E33.lat, aae_MY24[:,:], colors='black', levels=energylim)
py.clabel(fig, fontsize=9,inline=1, fontcolor='black', fmt='%1.0f')
ax4.set_xlabel('Latitude')
ax4.set_ylabel('Longitude')
fig2 = py.imread(datafile)
ax4.imshow(fig2, zorder=0, extent=[0, 360, 0, 60])
img=py.plot(longitude, latitude, 'ok',label='ROIs', color='blue',markersize=8)
img=py.plot(lon_crater, lat_crater, 'ok', color='yellow',markersize=8)
img=py.plot(lon_new, lat_new, 'ok', color='red', markersize=8)
ax4.grid(b=None, which='major', axis='both')

#fig=py.contourf(DDS_E33.lon, DDS_E33.lat, aae_MY24[:,:], cmap=py.cm.Greys, levels=energylim, alpha=.5)
#py.show()


# ADD NEW ROI based on wind energy 14 kW limit
ax4.text(234, 82, 'Hartwick+2021',color='tab:red',size=16)
ax4.add_patch(Rectangle((125,7),35,32,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((0,-7),77,57,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((33,-63),75,50,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((188,58),130,22,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((200,-25),66,80,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((275,-5),20,20,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((273,40),20,13,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((288,-60),29,22,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((168,-18),29,25,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((84,-4),27,25,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((136,-15),27,18,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((248,-49),27,19,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((320,-18),10,10,linewidth=4,edgecolor='red',facecolor='none'))
ax4.add_patch(Rectangle((316,1),10,10,linewidth=4,edgecolor='red',facecolor='none'))


# add landing sites
ax4.add_patch(Rectangle((99,-45),17,12, linewidth=2,edgecolor='tab:blue',facecolor='none')) #Holt+2008
#ax4.text(99-2, -45+14, 'Holt+2008',color='tab:blue', size=14)
ax4.add_patch(Rectangle((168,-35),9,10, linewidth=2,edgecolor='tab:orange',facecolor='none')) #Adeli+2018
#ax4.text(168-10, -35+12, 'Adeli+2018',color='tab:orange', size=14)
ax4.add_patch(Rectangle((300,-41),30,11, linewidth=2,edgecolor='tab:purple',facecolor='none')) #Berman+2018
#ax4.text(300, -41+13, 'Berman+2008',color='tab:purple', size=14)

ax4.add_patch(Rectangle((2,0),29,16, linewidth=2,edgecolor='tab:cyan',facecolor='none')) #Plaut+2009,Peterson+2018
#ax4.text(2, -7, 'Plaut+2009,Peterson+2018',color='tab:cyan', size=14)
ax4.add_patch(Rectangle((65,35),20,13, linewidth=2,edgecolor='tab:olive',facecolor='none')) #Stuurman+2016
#ax4.text(65-5, 49, 'Stuurman+2016',color='tab:olive', size=14)
#ax4.text(180, 27, 'Bramson+2015',color='tab:purple', size=14)

py.savefig(f'{PATH}/windenergy_fig{ct}.eps',dpi=200)

