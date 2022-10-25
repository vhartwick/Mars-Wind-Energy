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
bins = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13,13.5,14.,14.5,15,
        15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27.,27.5]
# Import Base Datasets
dataDIR = /lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/Data/FigureS13.nc'
DDS = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')

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
fig1=ax1.plot(DDS.areo[:-1],DDS.temp.sel(zagl=50,time_of_day_24=3.5)[:-1],color='tab:red',linestyle='dashed')
fig1=ax1.plot(DDS.areo[:-1],DDS.temp.sel(zagl=50,time_of_day_24=15.5)[:-1],color='tab:red')
axt = ax1.twinx()
axt.plot(DDS.areo[:-1],DDS.taudust_VIS.sel(lat=dm2lat,lon=dm2lon,time_of_day_24=3.5,method='nearest')[:-1],linestyle='dashed',color='tab:blue')
axt.plot(DDS.areo[:-1],DDS.taudust_VIS.sel(lat=dm2lat,lon=dm2lon,time_of_day_24=15.5,method='nearest')[:-1],color='tab:blue')
ax1.set_xlabel("$L_{s}$")
ax1.set_ylabel('[K]')
axt.set_ylabel(r'$[\tau]$')
ax1.text(0.95, 0.95, '(a)', ha='right', va='top',transform=ax1.transAxes, fontsize=12)
ax1.set_xlim([0,360])
ax1.text(50,232,'50m \n Temperature',color='tab:red')
ax1.text(265,240, 'Visible Dust \n Optical Depth',color='tab:blue')

# panel b - annual average vertical profile of rho/wind speed
ax2.xaxis.set_major_locator(MultipleLocator(0.001))
fig2  = ax2.plot(DDS.rho.mean('time').sel(time_of_day_24=3.5),DDS.zagl,marker='o',linestyle='dashed',color='tab:green')
ax2.plot(DDS.rho.mean('time').sel(time_of_day_24=15.5),DDS.zagl,marker='o',color='tab:green')
ax2.set_xlabel('$[kg/m^{3}$]')
axt = ax2.twiny()
axt.plot(DDS.temp.mean('time').sel(time_of_day_24=3.5), DDS.zagl,linestyle='dashed',marker='^',color='tab:red')
axt.plot(DDS.temp.mean('time').sel(time_of_day_24=15.5),DDS.zagl,marker='^',color='tab:red')
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
#fig3 = ax3.plot(DDS.uls90.mean('time'),DDS.vls90.mean('time'),marker='o')
fig3=ax3.scatter(DDS.uls90.mean('time'),DDS.vls90.mean('time'), c=DDS.time_of_day_24,cmap=py.cm.plasma)
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
ax.bar(DDS.wdls90.flatten(),DDS.wsls90.flatten(),normed=True, opening=0.8, edgecolor='white',bins=wsbins)

# panel e - Ls90 wind speed distribution
fig5 = ax5.bar(bins, DDS.percent_time_ls90*100, color='blue',edgecolor='black')
ax5.set_xlabel('Wind Speed [m/s]')
ax5.set_ylabel('Frequency [%]')
ax5.text(0.95, 0.95, '(e)', ha='right', va='top',transform=ax5.transAxes, fontsize=12)

# Ls=90 Label
py.gcf().text(0.05, 0.155, '$L_{s}=270\degree$', rotation='vertical')


# panel f - Ls270 u vs v
fig6=ax6.scatter(DDS.uls270.mean('time'),DDS.vls270.mean('time'),c=DDS.time_of_day_24, cmap=py.cm.plasma)
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
ax.bar(DDS.wdls270.flatten(),DDS.wsls270.flatten(),normed=True, opening=0.8, edgecolor='white',bins=wsbins)
ax.set_legend(bbox_to_anchor=(-0.5,-0.1),prop={'size':4})

# panel h - Ls270 wind speed distribution
fig8 = ax8.bar(bins, DDS.percent_time_ls270*100, color='blue',edgecolor='black')
ax8.set_xlabel('Wind Speed [m/s]')
ax8.set_ylabel('Frequency [%]')
ax8.text(0.95, 0.95, '(h)', ha='right', va='top',transform=ax8.transAxes, fontsize=12)

py.savefig(f'{PATH}/WindEnergy_HighRes_FigS{ct}.eps',dpi=300)
