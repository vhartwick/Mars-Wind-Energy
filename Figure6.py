# Figure 6: Python script for base analysis and figures for Hartwick+2022
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
import matplotlib.patches as mpatches

# Get Path
PATH = os.getcwd()

# Import Base Datasets
dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_solarwindpercent2500.nc'
DDS_PCT = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_solarwindpercent2500_season.nc'
DDS_PCT2 = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/lou/la2/vhartwic/FV3/verona/fv3_2.0_windenergy/fms_mars_MY24_highres/history/MY24_highres_savefile.nc'
DDS_swflx = xr.open_dataset(dataDIR, decode_times=False)

# define locations (lat/lon)
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

# -------------------------------------------------------------------------------
# FIGURE 6
# Fig. 6. We identify 13 new potential regions of interest where wind power
# exceeds 24kW for more than 40% of the Mars year (orange contours) or season
# (Ls90=light orange, Ls270=purple). Grey contours show regions where the
# combined solar and wind power exceeds 24kW more than 65% of the Mars year.
# Black and white boxes show newly identified ROI. Yellow circles show ROI that
# have been specifically flagged for future work based on our analysis (Table S2).
# Blue circles show crater ROI that may have significant topographic induced winds
# and should be studied at higher resolution. The background image shows results
# from the NASA Subsurface Water Ice Mapping (SWIM) 2.0 Global Products data
# release (8). Blue regions indicate ice with 1-5m below the surface.
# -------------------------------------------------------------------------------
ct = 6
print('Making Figure 6')

# find number of timesteps in file
time = len(DDS_swflx.time)*len(DDS_swflx.time_of_day_24)
levels = np.arange(0,100,5)
level_lim = [0,50,1000]


fig, axs = py.subplots(1,3, sharey=True, sharex=True, figsize=(12,5)) # define number of rows, columns
fig.subplots_adjust(wspace=0.05, hspace=0.01)

# panel a - wind
cs = axs[0].contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_wind24*100, cmap=py.cm.viridis,levels=levels)
cs2 = axs[0].contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_wind24*100, colors='k', levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[0].xaxis.set_major_locator(MultipleLocator(60))
axs[0].yaxis.set_major_locator(MultipleLocator(30))
axs[0].set_title('(a) Wind')
axs[0].set_xlabel('Longitude')
axs[0].set_ylabel('Latitude')

# panel b - solar
cs = axs[1].contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_solar2500*100, cmap=py.cm.viridis,levels=levels)
cs2 = axs[1].contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_solar2500*100,  colors='k',levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[1].set_title('(b) Solar, $\ 2500m^{2}$')
axs[1].set_xlabel('Longitude')

# panel c - total
cs = axs[2].contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_total2500*100, cmap=py.cm.viridis, levels=levels)
cs2 = axs[2].contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_total2500*100, colors='k',levels=levels[::3])
py.clabel(cs2, fontsize=9,inline=1, fmt='%1.0f')
axs[2].set_title('(c) Solar + Wind')
axs[2].set_xlabel('Longitude')

norm= matplotlib.colors.Normalize(vmin=cs.cvalues.min(), vmax=cs.cvalues.max())
sm = py.cm.ScalarMappable(norm=norm, cmap = cs.cmap)
sm.set_array([])
cbar_ax = py.gcf().add_axes([0.92, 0.125, 0.01, 0.755]) #left, bottom, width, height
clb = fig.colorbar(sm, cax=cbar_ax, orientation='vertical', label='% year')

py.savefig(f'{PATH}/WindEnergy_Highres_Fig{ct}.eps',dpi=300)

# panel d

levels=[39.9,40,100]

level_lim = [64.9,65,1000]
level_lim2 = [0,50,100]
fig, ax4 = py.subplots(1,1, sharey=True, figsize=(18,12)) # define number of rows, columns
fig.subplots_adjust(wspace=0.1, hspace=0.1)

fig1 = ax4.contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_total2500*100, levels=level_lim, colors='#bbbbbb')
fig1 = ax4.contour(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_total2500*100, levels=level_lim, colors='k')

fig2 = ax4.contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT2.percent_wind24_ls270*100, levels=levels, colors='#8a226a')
fig2 = ax4.contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT2.percent_wind24_ls90*100, levels=levels, colors='#f98e09')
fig2 = ax4.contourf(DDS_PCT.lon, DDS_PCT.lat, DDS_PCT.percent_wind24*100, levels=levels, colors='#e45a31')

ax4.xaxis.set_major_locator(MultipleLocator(60))
ax4.yaxis.set_major_locator(MultipleLocator(30))
ax4.set_xlim([0,360])
ax4.set_ylim([-90,90])
ax4.set_xlabel('Longitude', fontsize=14)
ax4.set_ylabel('Latitude', fontsize=14)

datafile = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/SWIM1-5_0180.png'
datafile2 = '/lou/la2/vhartwic/analysis/papers/windenergy/Hartwick2022/SWIM1-5_180360.png'
fig2 = py.imread(datafile)
fig3 = py.imread(datafile2)
ax4.imshow(fig2, zorder=0, extent=[180,360, -60, 60])
ax4.imshow(fig3, zorder=0, extent=[0,180,-60,60])

img=py.plot(longitude, latitude, 'ok',label='ROIs', markersize=12)
img=py.plot(lon_crater, lat_crater, 'o',color='#02d8e9',markersize=12, markeredgecolor='k',markeredgewidth=1.5, label='Crater ROI')
img=py.plot(lon_new, lat_new, 'o', color='#c9ff27', markersize=12, markeredgecolor='k',markeredgewidth=1.5, label='New ROI')
ax4.grid(visible=None, which='major', axis='both')

ax4.add_patch(Rectangle((22,-60),87,23,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((70,-37),32,19,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((167,-21),38,33,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((170,12),27,23,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((255,30),50,30,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((255,-63),60,25,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((308,-27),30,33,linewidth=4,edgecolor='w',facecolor='none'))
ax4.add_patch(Rectangle((265,-15),40,45,linewidth=4,edgecolor='w',facecolor='none'))

ax4.add_patch(Rectangle((5,60),120,23,linewidth=4,edgecolor='k',facecolor='none'))
ax4.add_patch(Rectangle((147,60),185,23,linewidth=4,edgecolor='k',facecolor='none'))
ax4.add_patch(Rectangle((205,-23),60,80,linewidth=4,edgecolor='white',facecolor='none'))


# First Legend
first_legend=ax4.legend(ncol=3,loc='lower center')
first_legend =ax4.legend(ncol=3, bbox_to_anchor=(0.5,0.05),loc='lower center')
ax4.add_artist(first_legend)

# second lengend
#p1 = mpatches.Rectangle((0,0),1,1,fc="#f98e09",alpha=0.5)
#handles = [p1,img]
#labels=['>40% of season, Ls=90', 'ROI','Crater ROI', 'New ROI']
#ax4.legend(handles,labels,ncol=4, loc='lower center')

p1 = mpatches.Patch(color='#bbbbbb',label='>65% of year, Wind+Solar')
p2 = mpatches.Patch(color='#e45a31', label='>40% of year, Wind')
p3 = mpatches.Patch(color='#8a226a',label='>40% of Ls270, Wind')
p4 = mpatches.Patch(color='#f98e09', label='>40% of Ls90, Wind')
ax4.legend(handles=[p1,p2,p3,p4],ncol=4,loc='lower center')
#proxy = [py.Rectagle(0,0),1,1,set_fc=['k','#bbbbbb','#8a226a','#f98e09']]
#ax4.legend(proxy,['>65% of year, Solar+Wind', '>40% of year, Wind', '>40% of season, Ls=270','>40% of season, Ls=90'])
#ax4.legend(nrow=2,ncol=3, loc='lower center')

py.savefig(f'{PATH}/WindEnergy_Highres_Fig{ct}b.pdf',dpi=300)

# Send out Data
var1 = xr.DataArray(DDS_PCT.percent_wind24, coords={'lon':DDS_PCT.lon, 'lat':DDS_PCT.lat},dims=['lat','lon'])
var1.name = 'percent_wind24'

var2 = xr.DataArray(DDS_PCT.percent_solar2500, coords={'lon':DDS_PCT.lon, 'lat':DDS_PCT.lat},dims=['lat','lon'])
var2.name = 'percent_solar2500'

var3 = xr.DataArray(DDS_PCT.percent_total2500, coords={'lon':DDS_PCT.lon, 'lat':DDS_PCT.lat},dims=['lat','lon'])
var3.name = 'percent_total2500'

var4 = xr.DataArray(DDS_PCT2.percent_wind24_ls270, coords={'lon':DDS_PCT.lon, 'lat':DDS_PCT.lat},dims=['lat','lon'])
var4.name = 'percent_wind24_ls270'

var5 = xr.DataArray(DDS_PCT2.percent_wind24_ls90, coords={'lon':DDS_PCT.lon, 'lat':DDS_PCT.lat},dims=['lat','lon'])
var5.name = 'percent_wind24_ls90'

NEW_DF = xr.merge([var1,var2,var3,var4,var5])
NEW_DF.to_netcdf(f'{PATH}/Data/Figure{ct}.nc')
