# python routine to plot environment at DM2 (35N,23E)
# author: victoria l. hartwick
# last modified: 6/29/22


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
from windrose import WindroseAxes

# Get Path
PATH = os.getcwd()

# Specify DM2 Lat/Lon
dm2lat = 35
dm2lon = 23

# Import Base Datasets
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/*.atmos_diurn_T_zagl_min.nc'
DDS_MY24 = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')

dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/*.atmos_diurn_T_temp_zagl.nc'
DDS = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')
temp = DDS.temp

print('Calculating Wind Direction')
u = DDS_MY24.ucomp.sel(zagl=50).sel(lat=35,lon=23,method='nearest').mean('time_of_day_24')
v = DDS_MY24.vcomp.sel(zagl=50).sel(lat=35,lon=24,method='nearest').mean('time_of_day_24')
ws = ((u**2) + (v**2))**0.5
ws = (u**2 + v**2)**0.5
wd = np.arctan(v/u)*(180/np.pi)
wd = np.where(wd > 0, wd, wd+360)

# note number of years
yrs = [0,1,2,3]

# locate cardinal seasons in Ls array
print('Finding Seasons')
areo = DDS_MY24.areo[:,0]
areo = np.mod(areo,360)

lsrange, lsrange1, lsrange2, lsrange3 = [], [], [], []
for i in yrs:
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

u = DDS_MY24.ucomp.sel(zagl=50).sel(lat=35,method='nearest').sel(lon=23,method='nearest')
v = DDS_MY24.vcomp.sel(zagl=50).sel(lat=35,method='nearest').sel(lon=23,method='nearest')

uls90,uls270 = u[ls90,:],u[ls270,:]
vls90,vls270 = v[ls90,:],v[ls270,:]

wsls90,wsls270 = (uls90**2 + vls90**2)**0.5, (uls270**2+vls270**2)**0.5
wdls90,wdls270 = np.arctan(vls90/uls90)*(180/np.pi), np.arctan(vls270/uls270)*(180/np.pi)
wdls90 = np.where(wdls90 > 0, wdls90, wdls90+360)
wdls270 = np.where(wdls270 > 0, wdls270, wdls270+360)
wsls90, wsls270 = np.array(wsls90),np.array(wsls270)
wdls90, wdls270 = np.array(wdls90), np.array(wdls270)

# other selected variables
rho = DDS.rho.sel(lat=dm2lat,lon=dm2lon,method='nearest')
temp = DDS.temp.sel(lat=dm2lat,lon=dm2lon,method='nearest')
dm2u = DDS_MY24.ucomp.sel(lat=dm2lat,lon=dm2lon,method='nearest')
dm2v = DDS_MY24.vcomp.sel(lat=dm2lat,lon=dm2lon,method='nearest')

# calculate seasonal wind speed distribution
print('Calculating Seasonal Wind Speed Distribution')
count_ls90, count_ls270 = [],[]
bins = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13,13.5,14.,14.5,15,
        15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27.,27.5]
bins_top = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13,13.5,14.,14.5,15,
        15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27.,27.5,28]

for x,y in zip(bins, bins_top):   # loop through wind speed bins
	tmpls90 = np.count_nonzero(np.where(np.logical_and(wsls90 >=x,wsls90 < y)))
	tmpls270 = np.count_nonzero(np.where(np.logical_and(wsls270 >=x,wsls270 <y)))
	count_ls90 += [tmpls90]  # number of original velocity measurements in each bin, at atm_lev
	count_ls270 += [tmpls270]  # number of original velocity measurements in each bin, at atm_lev

# calculate percent of time operating in each interval
percent_time_ls90 = np.divide(count_ls90, np.sum(count_ls90))
percent_time_ls270 = np.divide(count_ls270, np.sum(count_ls270))

# Plot and Save Figure
ct = 13
print('Making Figure S13')
fig, axs = py.subplots(4,4, sharey=True, figsize=(12,10)) # define number of rows, columns
gs = gridspec.GridSpec(4,4, wspace=0.4,hspace=0.4) # define grid matching number of rows, columns

ax1 = py.subplot(gs[0:2,0:3]) # define size of each panel row x column
ax2 = py.subplot(gs[0:2,3])
ax3 = py.subplot(gs[2,0])
ax4 = py.subplot(gs[2,1])
ax5 = py.subplot(gs[2,2:])
ax6 = py.subplot(gs[3,0])
ax7 = py.subplot(gs[3,1])
ax8 = py.subplot(gs[3,2:])


# panel a - 3am/3pm surface temperature, dust optical depth
fig1=ax1.plot(areo[:-1],temp.sel(zagl=50,time_of_day_24=3.5)[:-1],color='tab:red',linestyle='dashed')
fig1=ax1.plot(areo[:-1],temp.sel(zagl=50,time_of_day_24=15.5)[:-1],color='tab:red')
axt = ax1.twinx()
axt.plot(areo[:-1],DDS_MY24.taudust_VIS.sel(lat=dm2lat,lon=dm2lon,time_of_day_24=3.5,method='nearest')[:-1],linestyle='dashed',color='tab:blue')
axt.plot(areo[:-1],DDS_MY24.taudust_VIS.sel(lat=dm2lat,lon=dm2lon,time_of_day_24=15.5,method='nearest')[:-1],color='tab:blue')
ax1.set_xlabel("$L_{s}$")
ax1.set_ylabel('[K]')
axt.set_ylabel(r'$[\tau]$')
ax1.text(0.95, 0.95, '(a)', ha='right', va='top',transform=ax1.transAxes, fontsize=12)
ax1.set_xlim([0,360])
ax1.text(50,232,'50m \n Temperature',color='tab:red')
ax1.text(265,240, 'Visible Dust \n Optical Depth',color='tab:blue')

# panel b - annual average vertical profile of rho/wind speed
ax2.xaxis.set_major_locator(MultipleLocator(0.001))
fig2  = ax2.plot(rho.mean('time').sel(time_of_day_24=3.5),DDS.zagl,marker='o',linestyle='dashed',color='tab:green')
ax2.plot(rho.mean('time').sel(time_of_day_24=15.5),DDS.zagl,marker='o',color='tab:green')
ax2.set_xlabel('$[kg/m^{3}$]')
axt = ax2.twiny()
axt.plot(temp.mean('time').sel(time_of_day_24=3.5), DDS.zagl,linestyle='dashed',marker='^',color='tab:red')
axt.plot(temp.mean('time').sel(time_of_day_24=15.5),DDS.zagl,marker='^',color='tab:red')
axt.set_ylabel('Altitude [km]')
axt.set_xlabel ('[K]')
ax2.text(0.95, 0.95, '(b)', ha='right', va='top',transform=ax2.transAxes, fontsize=12)
ax2.set_ylim([5,100])
ax2.set_xlim([0.007, 0.009])
axt.set_xlim([180,250])
ax2.text(20,80,r"$\rho$",color='tab:blue')
ax2.text(26,27,'wind speed',color='tab:green')
ax2.yaxis.tick_right()

# Ls=90 Label
py.gcf().text(0.05, 0.37, '$L_{s}=90\degree$', rotation='vertical')

# panel c - Ls90 u vs v
#fig3 = ax3.plot(uls90.mean('time'),vls90.mean('time'),marker='o')
fig3=ax3.scatter(uls90.mean('time'),vls90.mean('time'), c=DDS_MY24.time_of_day_24,cmap=py.cm.plasma)
ax3.set_xlabel('U [m/s]')
ax3.set_ylabel('V [m/s]')
ax3.text(0.95, 0.95, '(c)', ha='right', va='top',transform=ax3.transAxes, fontsize=12)

# panel d - Ls90 wind rose
ax4.axis('off')
wsbins=np.arange(0,30,6)
gp = gs.get_grid_positions(fig)  # [bottom, top, left, right]
rect = [gp[2][1], gp[2][1],
                gp[3][0]-gp[2][0],
                gp[1][0]-gp[0][0]]  # [left, bottom, width, height]
ax = WindroseAxes(fig, rect)
fig.add_axes(ax)
ax.bar(wdls90.flatten(),wsls90.flatten(),normed=True, opening=0.8, edgecolor='white',bins=wsbins)

# panel e - Ls90 wind speed distribution
fig5 = ax5.bar(bins, percent_time_ls90*100, color='blue',edgecolor='black')
ax5.set_xlabel('Wind Speed [m/s]')
ax5.set_ylabel('Frequency [%]')
ax5.text(0.95, 0.95, '(e)', ha='right', va='top',transform=ax5.transAxes, fontsize=12)

# Ls=90 Label
py.gcf().text(0.05, 0.155, '$L_{s}=270\degree$', rotation='vertical')


# panel f - Ls270 u vs v
fig6=ax6.scatter(uls270.mean('time'),vls270.mean('time'),c=DDS_MY24.time_of_day_24, cmap=py.cm.plasma)
ax6.set_xlabel('U [m/s]')
ax6.set_ylabel('V [m/s]')
ax6.text(0.95, 0.95, '(f)', ha='right', va='top',transform=ax6.transAxes, fontsize=12)

# panel g - Ls270 wind rose
ax7.axis('off')
gp = gs.get_grid_positions(fig)  # [bottom, top, left, right]
rect = [gp[2][1], gp[0][3],
                gp[3][0]-gp[2][0],
                gp[1][0]-gp[0][0]]  # [left, bottom, width, height]
ax = WindroseAxes(fig, rect)
fig.add_axes(ax)
ax.bar(wdls270.flatten(),wsls270.flatten(),normed=True, opening=0.8, edgecolor='white',bins=wsbins)
ax.set_legend(bbox_to_anchor=(-0.5,-0.1),prop={'size':4})

# panel h - Ls270 wind speed distribution
fig8 = ax8.bar(bins, percent_time_ls270*100, color='blue',edgecolor='black')
ax8.set_xlabel('Wind Speed [m/s]')
ax8.set_ylabel('Frequency [%]')
ax8.text(0.95, 0.95, '(h)', ha='right', va='top',transform=ax8.transAxes, fontsize=12)

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)

# Save Data
var1 = xr.DataArray(areo, coords={'time':DDS_MY24.time},dims=['time'])
var2 = xr.DataArray(temp, coords={'time':DDS_MY24.time,'zagl':DDS_MY24.zagl,'time_of_day_24':DDS_MY24.time_of_day24},dims=['time','time_of_day_24','zagl'])
var3 = xr.DataArray(rho, coords={'time':DDS_MY24.time,'zagl':DDS_MY24.zagl,'time_of_day_24':DDS_MY24.time_of_day24},dims=['time','time_of_day_24','zagl'])
NEW_DF = xr.Dataset({'areo':var1,'temp':var2,'taudust_VIS':DDS_MY24.taudust_VIS,'rho':var3,
	'uls90':uls90,'vls90':vls90,'wdls90':wdls90,'wsls90':wsls90,'percent_time_ls90':percent_time_ls90,
	'uls270':uls270,'vls270':vls270,'wdLs270':wdLs270,'wsLs270':wsLs270,'percent_time_ls270':percent_time_ls270})

NEW_DF.to_netcdf(f'{PATH}/Data/FigureS{ct}.nc')
