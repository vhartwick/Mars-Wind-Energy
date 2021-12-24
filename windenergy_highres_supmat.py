
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

# decide which figures to plot
i = 2

# read in first dataset
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/*.atmos_diurn_T_zagl_min.nc'
DDS_SY1 = xr.open_mfdataset(dataDIR, decode_times=False,combine='by_coords')
   
dataDIR = '/Users/vhartwic/Documents/codes/VHPlots/Data/fv3_2.0_windenergy/MY24_highres/02094.fixed.nc'
DS_fixed = xr.open_dataset(dataDIR, decode_times=False)

dataDIR = '/Users/vhartwic/Downloads/Enercon_330kW_power_MY24_HIGHRES.nc'
DDS_E33 = xr.open_dataset(dataDIR, decode_times=False)

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
#taudust_annual_ave_SY1 = DAS_SY1.taudust_VIS.mean('time')

wpd_time_ave_SY1 = wpd_SY1.mean('time_of_day_24')
wpd_diurn_ave_SY1 = wpd_SY1.mean('time')
swflx_diurn_ave_SY1 = swflx_diurn_SY1.mean('time')

# calculate solar power
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

#---------------------------------------------------------------------------------
# FIGURE S1
#---------------------------------------------------------------------------------
if i == 1:
	print('Making Figure S1')
	ct = 1

	# plot specific location
	tmp_MY24 = wpd_SY1.mean('time_of_day_24')
	tmp2_MY24 = DDS_SY1.swflx.mean('time_of_day_24')
	tmp2_MY24 = tmp2_MY24.mean('lon')
	areo_MY24 = DDS_SY1.areo[:,0]-360*3
	fig, ax = py.subplots(5,1, sharex=True, figsize=(12,10))
	fig.subplots_adjust(wspace=0.05)

	fig1 = ax[0].plot(areo_MY24, tmp2_MY24.sel(lat=70,  method='nearest'), linewidth=3,linestyle='--')
	fig1 = ax[0].plot(areo_MY24, tmp_MY24.sel(lat=70, zagl=50, method='nearest'))
	ax[0].set_ylabel('$\ [W/m^{2}]$')
	ax[0].xaxis.set_major_locator(MultipleLocator(60))
	ax[0].set_xlim([0,360])

	fig1 = ax[1].plot(areo_MY24, tmp2_MY24.sel(lat=50,  method='nearest'),linewidth=3,linestyle='--')
	fig1 = ax[1].plot(areo_MY24, tmp_MY24.sel(lat=50, zagl=50, method='nearest'))
	ax[1].set_ylabel('$\ [W/m^{2}]$')

	fig1 = ax[2].plot(areo_MY24, tmp2_MY24.sel(lat=35,  method='nearest'),linewidth=3,linestyle='--')
	fig1 = ax[2].plot(areo_MY24, tmp_MY24.sel(lat=35, zagl=50, method='nearest'))

	#ax[2].set_title('35N, all lon')
	ax[2].set_ylabel('$\ [W/m^{2}]$')


	fig2 = ax[3].plot(areo_MY24, tmp2_MY24.sel(lat=-43,  method='nearest'),linewidth=3,linestyle='--')
	fig2 = ax[3].plot(areo_MY24, tmp_MY24.sel(lat=-43, zagl=50, method='nearest'))

	#ax[3].set_title('43S, all lon, Holt 2008')
	ax[3].set_ylabel('$\ [W/m^{2}]$')

	#a = np.mean(tmp_SY1[:,:,:,24:33],axis=3)   # average between lon = 90-120 E
	fig2 = ax[4].plot(areo_MY24, tmp2_MY24.sel(lat=-50,  method='nearest'),linewidth=3,linestyle='--')
	fig2 = ax[4].plot(areo_MY24, tmp_MY24.sel(lat=-50, zagl=50, method='nearest'))

	#ax[4].set_title('50S, all lon, Holt 2008')
	ax[4].set_ylabel('$\ [W/m^{2}]$')
	ax[4].set_xlabel('$\ L_{s}$')

	py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}.epg', dpi=200)

#---------------------------------------------------------------------------------
# FIGURE S2
#---------------------------------------------------------------------------------
print('Making FigureS2')
ct=2

levels = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]
levels = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
#levels = [200,250,300,350,400,450]
levels_swflx = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]

fig, axs = py.subplots(4,6, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)
    
fig.suptitle('Annual Average Wind Profile Density $\ [W/m^{2}$], 50m')
fig.subplots_adjust(hspace=0.35)
axs = axs.ravel()
for i in range(24):

    fig1 = axs[i].contourf(DDS_SY1.lon,DDS_SY1.lat,wpd_diurn_ave_SY1.sel(zagl=50,time_of_day_24=i+0.5),levels=levels,cmap=py.cm.viridis)
    axs[i].set_title(str(i)+'LT')

    if i >=18:
        axs[i].set_xlabel('Longitude')
        
    #add colorabar
    #py.colorbar(fig1,ax=axs[i],orientation='horizontal')
    if i == 0 or i ==6 or i==12 or i==18:
        axs[i].set_ylabel('Latitude')


py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}.eps', dpi=200)
# also look at solar power
fig, axs = py.subplots(4,6, sharey=True, sharex=True, figsize=(12,5))
fig.subplots_adjust(wspace=0.05)


    
fig.suptitle('Annual Average Solar Power $\ [W/m^{2}$]')

for i in range(24):

    fig1 = axs[i].contourf(DDS_MYSY1.lon,DDS_MYSY1.lat,swflx_diurn_ave_MY24.sel(time_of_day_24=i+0.5),levels=levels,cmap=py.cm.viridis)
    axs[i].set_title(str(i)+'LT')
    if i>=18:
        axs[i].set_xlabel('Longitude')
        
    #add colorabar
    if i == 0 or i==6 or i==12 or i ==18:
        axs[i].set_ylabel('Latitude')

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}b.eps', dpi=200)

#---------------------------------------------------------------------------------
# FIGURE S3
#---------------------------------------------------------------------------------
print('Making Figure S3')
ct =3
# MAKE SURE TO CHECK EFFICIENCY FACTOR FOR GENERAL TURBINE AND CHANGE LABELS 

###### goal here is to count the number of times velocities at a particular location fall in each bin
# as a fraction of the total time. This is pretty easy for the average files but artificially drops down
# wind speed because of diurnal average. I can use the diurnal wind speed and first find the counts for each local time
# and then do soemthing else
# calculate based on power (preferred)

# this would be to count specifically in intervals
count_SY1, count_SY2, count2_SY1, count2_SY2 = [], [], [], []
bins = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13,13.5,14.,14.5,15,
        15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27.,27.5]
bins_top = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8.,8.5,9.,9.5,10.,10.5,11.,11.5,12.,12.5,13,13.5,14.,14.5,15,
        15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27.,27.5,28]

    
for x,y in zip(bins, bins_top):   # loop through wind speed bins
    tmp = np.count_nonzero(np.where(np.logical_and(DDS_E33.veladj_diurn_SY1.sel(zagl=atm_lev, method='nearest') >=x, 
                                                   DDS_E33.veladj_diurn_SY1.sel(zagl=atm_lev, method='nearest')  < y)))
    tmp2 = np.count_nonzero(np.where(np.logical_and(vel_diurn_SY1.sel(zagl=atm_lev, method='nearest') >=x, 
                                                    vel_diurn_SY1.sel(zagl=atm_lev, method='nearest') < y)))
   
    count_SY1 += [tmp]    # number of adjusted velocity measurements in each bin, at atm_lev
    count2_SY1 += [tmp2]  # number of original velocity measurements in each bin, at atm_lev


# calculate percent of time operating in each interval
percent_time_SY1 = np.divide(count_SY1,np.sum(count_SY1))
percent_time2_SY1 = np.divide(count2_SY1, np.sum(count2_SY1))


# plot
fig, ax1 = py.subplots(1,1, figsize=(12,10))
fig.subplots_adjust(wspace=0.05)

fig1 = ax1.bar(bins, percent_time2_SY1*100, color='white', edgecolor='black')
fig1 = ax1.bar(bins,percent_time_SY1*100, hatch='/', color='grey', edgecolor='black')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

fig1 = ax2.plot(bins,power333kwopt, linewidth=3, label='Enercon E33, 330kW')
ax2.set_ylim(0,750)
ax2.plot(bins,genE333_powercurve*1e-3, linewidth=3, color='purple', label='30% Efficiency, 330 kW')
ax2.plot(bins,power20kwopt*10, linewidth=3, label='Jacobs 31-20, 20kW x 10')
ax2.plot(bins_5mw,power5mwopt*0.1, linewidth=3, label='NREL 5MW x 0.1')
ax2.set_title('Power Curve & Global Wind Speed Probability Density', fontsize=14)
ax2.set_ylabel('Turbine Power Output [kW]')
ax2.legend(loc='upper right')
   
ax1.set_ylabel('Wind Speed Distribution [%]')
ax1.set_xlabel('Wind Speed [m/s]')

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}.eps',dpi=200)

#---------------------------------------------------------------------------------
# FIGURE S4
#---------------------------------------------------------------------------------
print('Making Figure S4')
ct =4

# Calculate Load Duration Curves for Global Year & Season

# calculate based on power (preferred)
p333 = (0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
p333top = (3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)

# calculate global mean
zonal_ave_SY1 = np.mean(DDS_E33.pow_diurn_SY1, axis=3)

# make sure to weight by area
latr = np.deg2rad(DDS_E33.lat)
weights = np.cos(latr)
global_ave_SY1 = np.average(zonal_ave_SY1, axis=2, weights=weights)

# select seasons
global_ave_ls90_SY1 = global_ave_SY1[ls90,:]
global_ave_ls180_SY1 = global_ave_SY1[ls180,:]
global_ave_ls270_SY1 = global_ave_SY1[ls270,:]
global_ave_ls0_SY1 = global_ave_SY1[ls0,:]


power_range =(0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)



# what we really want is the amount of time capacity factor is x or greater so...
count_SY1, countls90_SY1, countls180_SY1, countls270_SY1, countls0_SY1 = [], [], [], [], []


for x in power_range:   # loop through 

    # SY1
    count_SY1 += [np.count_nonzero(np.where(global_ave_SY1[:,:] >= x))]
    countls90_SY1 += [np.count_nonzero(np.where(global_ave_ls90_SY1[:,:] >= x))]
    countls180_SY1 += [np.count_nonzero(np.where(global_ave_ls180_SY1[:,:] >= x))]
    countls270_SY1 += [np.count_nonzero(np.where(global_ave_ls270_SY1[:,:] >= x))]
    countls0_SY1 += [np.count_nonzero(np.where(global_ave_ls0_SY1[:,:] >= x))]
    
    
# calculate percent of time operating in each interval
percent_time_gave_SY1 = np.divide(count_SY1,count_SY1[0])*100
percent_time_ls90_SY1 = np.divide(countls90_SY1,countls90_SY1[0])*100
percent_time_ls180_SY1 = np.divide(countls180_SY1,countls180_SY1[0])*100
percent_time_ls270_SY1 = np.divide(countls270_SY1,countls270_SY1[0])*100
percent_time_ls0_SY1 = np.divide(countls0_SY1,countls0_SY1[0])*100


# calculate capacity factor (power interval / total)
percent_max_gave = np.divide(power_range, max(power_range))*100
#print(percent_max)

# plot
fig1, ax2 = py.subplots(constrained_layout=True)
CS = ax2.plot(percent_time_gave_SY1,percent_max_gave, linewidth=2, label='Annual Average')
ax2.plot(percent_time_ls90_SY1, percent_max_gave, linewidth=2, label='Ls90')
ax2.plot(percent_time_ls180_SY1, percent_max_gave, linewidth=2, label='Ls180')
ax2.plot(percent_time_ls270_SY1, percent_max_gave, linewidth= 2, label='Ls270')
ax2.plot(percent_time_ls0_SY1, percent_max_gave, linewidth=2, label='Ls0')
py.ylim(0,8)
py.xlim(0.,100)

ax2.add_patch(Rectangle((0,0),50,100, linewidth=1,edgecolor='tab:olive',facecolor='none', hatch='/')) #Stuurman+2016
ax2.set_title('Global Load Duration Curve')
ax2.set_xlabel('% of time')
ax2.set_ylabel('% capacity factor')
#py.plot(percent_time, percent_max)

#ax2.axhline(3.8, color='orange', linestyle='--')
py.legend()
print(percent_time_gave_SY1)
print(percent_max_gave)

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}a.eps',dpi=200)

### NEED TO SPECIFY INDEX FOR LAT/LON OF ROIS
# for now pick ROIs
#pow_5mw_diurn = np.where(np.logical_and(veladj_diurn >=x, veladj_diurn < y),z,pow_5mw_diurn)

lon_roi = [22.117202268431015, 22.797731568998103, 17.0132325141777,47.97731568998114]
lat_roi = [36.702439024390266, 34.75609756097563, 29.195121951219527,37.81463414634148]
roi_name = ['Deuteronilus Mensae', ' ', 'Ismenius Lacus', 'Protonilus Mensae']
color = ['red','purple','green']

# locate indeces of ROI lon/lat
lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_E33.lon-tmp[i])).argmin()
    lonlim_idx += [idx]
    
tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_E33.lat-tmp[i])).argmin()
    latlim_idx += [idx]

local_power_SY1 = DDS_E33.pow_diurn_SY1

# this plots the power versus time at that location
#py.plot(DS.areo, loc)

# calculate based on power (preferred)
p333 = (0.,3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)
p333top = (3., 5., 9., 14., 22., 30., 43., 55., 74., 92., 115., 138., 167., 196., 223., 250., 273., 293., 308., 320., 330., 335.)

# count number of times turbine is operating at or above power levels
count_roi1_SY1,count_roi2_SY1,count_roi3_SY1,count_roi4_SY1 = [],[],[],[]

for x in power_range:   # loop through 
    tmp_roi1_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_roi1_SY1 += [tmp_roi1_SY1]


    tmp_roi2_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[1],lonlim_idx[1]] >= x))
    count_roi2_SY1 +=[tmp_roi2_SY1]


    tmp_roi3_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[2],lonlim_idx[2]] >= x))
    count_roi3_SY1 +=[tmp_roi3_SY1]
    
    tmp_roi4_SY1 = np.count_nonzero(np.where(local_power_SY1[:,:,latlim_idx[3],lonlim_idx[3]] >= x))
    count_roi4_SY1 +=[tmp_roi4_SY1]

# calculate percent of time turbine is operating at each power interval
percent_time_roi1_SY1 = np.divide(count_roi1_SY1, count_roi1_SY1[0])*100
percent_time_roi2_SY1 = np.divide(count_roi2_SY1, count_roi2_SY1[0])*100
percent_time_roi3_SY1 = np.divide(count_roi3_SY1, count_roi3_SY1[0])*100
percent_time_roi4_SY1 = np.divide(count_roi4_SY1, count_roi4_SY1[0])*100


# calculate capacity factor (power interval / total)
percent_max = np.divide(power_range, max(power_range))*100
#print(percent_max)

# remove first integer
#percent_max=percent_max[1:]
#print(percent_time.shape)

# plot
fig1, ax2 = py.subplots(constrained_layout=True)


py.xlim(0,100)
py.ylim(0,50)
ax2.set_title('Annual Average Load Duration Curve')
ax2.set_xlabel('% of time')
ax2.set_ylabel('% capacity')

py.plot(percent_time_gave_SY1,percent_max_gave, linewidth=2, label='Annual Average')
py.plot(percent_time_roi1_SY1,percent_max, 'red',linewidth=2, label='Deuteronilus Mensae')
py.plot(percent_time_roi3_SY1, percent_max, 'purple',linewidth=2, label='Ismenius Lacus')
py.plot(percent_time_roi4_SY1, percent_max, 'green', linewidth=2, label='Protonilus Mensae')


ax2.add_patch(Rectangle((0,0),50,100, linewidth=1,edgecolor='tab:olive',facecolor='none', hatch='/')) #Stuurman+2016
py.legend()

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}b.eps',dpi=200)

# Finally Look at One Location for all Seasons
lon_roi = [47.97731568998114]
lat_roi = [37.81463414634148]
roi_name = ['Protonilus Mensae']

# locate indeces of ROI lon/lat
lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_E33.lon-tmp[i])).argmin()
    lonlim_idx += [idx]
    
tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_E33.lat-tmp[i])).argmin()
    latlim_idx += [idx]

local_power_SY1 = DDS_E33.pow_diurn_SY1


# this plots the power versus time at that location

# count number of times turbine is operating at or above power levels

count_ls90_SY1,count_ls180_SY1,count_ls270_SY1,count_ls0_SY1 = [],[],[],[]
for x in power_range:   # loop through 
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls90,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls90_SY1 += [tmp_SY1]
    
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls180,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls180_SY1 +=[tmp_SY1]

    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls270,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls270_SY1 +=[tmp_SY1]

    
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls0,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls0_SY1 +=[tmp_SY1]
    
# calculate percent of time turbine is operating at each power interval
ls90_SY1 = np.divide(count_ls90_SY1, count_ls90_SY1[0])*100
ls180_SY1 = np.divide(count_ls180_SY1, count_ls180_SY1[0])*100
ls270_SY1 = np.divide(count_ls270_SY1, count_ls270_SY1[0])*100
ls0_SY1 = np.divide(count_ls0_SY1, count_ls0_SY1[0])*100



# calculate capacity factor (power interval / total)
percent_max = np.divide(power_range, max(power_range))*100

# plot
fig1, ax2 = py.subplots(constrained_layout=True)

py.xlim(0,100)
py.ylim(0,50)
ax2.set_title('Ismenius Lacus Load Duration Curve [330kW]')
ax2.set_title('Protonilus Mensae Load Duration Curve')
ax2.set_xlabel('% of time')
ax2.set_ylabel('% capacity')
#ax2.set_yscale('log')
py.plot(percent_time_roi4_SY1, percent_max, label='Annual Average')
#py.plot(percent_time_gave_SY1,percent_max_gave, label='Annual Average')
py.plot(ls0_SY1, percent_max, 'green', label='Ls0')
py.plot(ls90_SY1, percent_max, 'orange', label='Ls90')
py.plot(ls180_SY1,percent_max, 'red', label='Ls180')
py.plot(ls270_SY1, percent_max, 'purple', label = 'Ls270')


py.legend()

#py.plot(percent_time_gave_SY2,percent_max_gave, linestyle='dotted')
#py.plot(ls180_SY2,percent_max, 'red', linestyle='dotted')
#py.plot(ls270_SY2, percent_max, 'purple', linestyle='dotted')
#py.plot(ls0_SY2, percent_max, 'green', linestyle='dotted')
#py.plot(ls90_SY2, percent_max, 'orange', linestyle='dotted')
ax2.add_patch(Rectangle((0,0),50,100, linewidth=1,edgecolor='tab:olive',facecolor='none', hatch='/')) #Stuurman+2016

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}c.eps',dpi=200)

# Finally Look at One Location for all Seasons
# ismenius lacus quadrangle 0-60E, 30-65 N
lon_roi = [17.0132325141777]
lat_roi = [29.195121951219527]
roi_name = ['Ismenius Lacus']

# locate indeces of ROI lon/lat
lonlim_idx, latlim_idx = [], []
tmp = np.array(lon_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_E33.lon-tmp[i])).argmin()
    lonlim_idx += [idx]
    
tmp = np.array(lat_roi)
for i in range(0,len(tmp)):
    idx = (np.abs(DDS_E33.lat-tmp[i])).argmin()
    latlim_idx += [idx]

local_power_SY1 = DDS_E33.pow_diurn_SY1


# this plots the power versus time at that location

# count number of times turbine is operating at or above power levels

count_ls90_SY1,count_ls180_SY1,count_ls270_SY1,count_ls0_SY1 = [],[],[],[]
for x in power_range:   # loop through 
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls90,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls90_SY1 += [tmp_SY1]
    
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls180,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls180_SY1 +=[tmp_SY1]

    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls270,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls270_SY1 +=[tmp_SY1]

    
    tmp_SY1 = np.count_nonzero(np.where(local_power_SY1[ls0,:,latlim_idx[0],lonlim_idx[0]] >= x))
    count_ls0_SY1 +=[tmp_SY1]
    
# calculate percent of time turbine is operating at each power interval
ls90_SY1 = np.divide(count_ls90_SY1, count_ls90_SY1[0])*100
ls180_SY1 = np.divide(count_ls180_SY1, count_ls180_SY1[0])*100
ls270_SY1 = np.divide(count_ls270_SY1, count_ls270_SY1[0])*100
ls0_SY1 = np.divide(count_ls0_SY1, count_ls0_SY1[0])*100



# calculate capacity factor (power interval / total)
percent_max = np.divide(power_range, max(power_range))*100

# plot
fig1, ax2 = py.subplots(constrained_layout=True)

py.xlim(0,100)
py.ylim(0,50)
ax2.set_title('Ismenius Lacus Load Duration Curve')

ax2.set_xlabel('% of time')
ax2.set_ylabel('% capacity')
#ax2.set_yscale('log')
py.plot(percent_time_roi4_SY1, percent_max, label='Annual Average')
#py.plot(percent_time_gave_SY1,percent_max_gave, label='Annual Average')
py.plot(ls0_SY1, percent_max, 'green', label='Ls0')
py.plot(ls90_SY1, percent_max, 'orange', label='Ls90')
py.plot(ls180_SY1,percent_max, 'red', label='Ls180')
py.plot(ls270_SY1, percent_max, 'purple', label = 'Ls270')


py.legend()

#py.plot(percent_time_gave_SY2,percent_max_gave, linestyle='dotted')
#py.plot(ls180_SY2,percent_max, 'red', linestyle='dotted')
#py.plot(ls270_SY2, percent_max, 'purple', linestyle='dotted')
#py.plot(ls0_SY2, percent_max, 'green', linestyle='dotted')
#py.plot(ls90_SY2, percent_max, 'orange', linestyle='dotted')
ax2.add_patch(Rectangle((0,0),50,100, linewidth=1,edgecolor='tab:olive',facecolor='none', hatch='/')) #Stuurman+2016

print(ls0_SY1)
print(percent_max)

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}d.eps',dpi=200)

#-----------------------------------------------------------------------------
# FIGURE S5
#-----------------------------------------------------------------------------
print('Making Figure S5')
# plot all locations with diurnal average power greater than some limit (ALL THE TIME)
# choose locations with energy (kW) always greater than the minimum energy requirements (24) 
# or the average greater than the mim energy requirement

# diurnal mean
tmp_SY1 = DDS_E33.pow_diurn_SY1.mean('time_of_day_24')


lat_sav24_ls90, lon_sav24_ls90, lat_sav24_ls270, lon_sav24_ls270= [], [], [], []
lat_sav15_ls90, lon_sav15_ls90, lat_sav15_ls270, lon_sav15_ls270= [], [], [], []
lat_sav2_ls90, lon_sav2_ls90, lat_sav2_ls270, lon_sav2_ls270= [], [], [], []


for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls90 = tmp_SY1[ls90,i,j]
        if np.any(tmp_ls90 < 24):
            continue
        else:
            lat_sav24_ls90+= [DDS_E33.lat[i]]
            lon_sav24_ls90 +=[DDS_E33.lon[j]]

            
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls270 = tmp_SY1[ls270,i,j]
        if np.any(tmp_ls270 < 24):
            continue
        else:
            lat_sav24_ls270+= [DDS_E33.lat[i]]
            lon_sav24_ls270 +=[DDS_E33.lon[j]]
            

# 15 kW limit
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls90 = tmp_SY1[ls90,i,j]
        if np.any(tmp_ls90 < 15):
            continue
        else:
            lat_sav15_ls90+= [DDS_E33.lat[i]]
            lon_sav15_ls90 +=[DDS_E33.lon[j]]
        
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls270 = tmp_SY1[ls270,i,j]
        if np.any(tmp_ls270 < 15):
            continue
        else:
            lat_sav15_ls270+= [DDS_E33.lat[i]]
            lon_sav15_ls270 +=[DDS_E33.lon[j]]

# 2 kW limit            
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls90 = tmp_SY1[ls90,i,j]
        if np.any(tmp_ls90 < 2):
            continue
        else:
            lat_sav2_ls90+= [DDS_E33.lat[i]]
            lon_sav2_ls90 +=[DDS_E33.lon[j]]
        
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls270 = tmp_SY1[ls270,i,j]
        if np.any(tmp_ls270 < 2):
            continue
        else:
            lat_sav2_ls270+= [DDS_E33.lat[i]]
            lon_sav2_ls270 +=[DDS_E33.lon[j]]
# genesate day and night masks
m = np.ma.masked_where(DDS_MYSY1.swflx_MY24>0,DDS_E33.pow_diurn_SY1)
print(m.shape)
E33_night = np.mean(m,axis=1)
print(E33_night.shape)

m = np.ma.masked_where(DDS_SY1.swflx_MY24==0,DDS_E33.pow_diurn_SY1)
E33_day = np.mean(m,axis=0)
E33_day = np.mean(E33_day,axis=0)

# plot all locations with diurnal average power greater than some limit (ALL THE TIME)
# choose locations with energy (kW) always greater than the minimum energy requirements (24) 
# or the average greater than the mim energy requirement

tmp_SY1 = E33_night


lat_sav24_ls90pm, lon_sav24_ls90pm, lat_sav24_ls270pm, lon_sav24_ls270pm= [], [], [], []
lat_sav15_ls90pm, lon_sav15_ls90pm, lat_sav15_ls270pm, lon_sav15_ls270pm= [], [], [], []
lat_sav2_ls90pm, lon_sav2_ls90pm, lat_sav2_ls270pm, lon_sav2_ls270pm= [], [], [], []


for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls90 = tmp_SY1[ls90,i,j]
        if np.any(tmp_ls90 < 24):
            continue
        else:
            lat_sav24_ls90pm+= [DDS_E33.lat[i]]
            lon_sav24_ls90pm +=[DDS_E33.lon[j]]

            
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls270 = tmp_SY1[ls270,i,j]
        if np.any(tmp_ls270 < 24):
            continue
        else:
            lat_sav24_ls270pm+= [DDS_E33.lat[i]]
            lon_sav24_ls270pm +=[DDS_E33.lon[j]]
            

# 15 kW limit
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls90 = tmp_SY1[ls90,i,j]
        if np.any(tmp_ls90 < 15):
            continue
        else:
            lat_sav15_ls90pm+= [DDS_E33.lat[i]]
            lon_sav15_ls90pm +=[DDS_E33.lon[j]]
        
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls270 = tmp_SY1[ls270,i,j]
        if np.any(tmp_ls270 < 15):
            continue
        else:
            lat_sav15_ls270pm+= [DDS_E33.lat[i]]
            lon_sav15_ls270pm +=[DDS_E33.lon[j]]

# 2 kW limit            
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls90 = tmp_SY1[ls90,i,j]
        if np.any(tmp_ls90 < 2):
            continue
        else:
            lat_sav2_ls90pm+= [DDS_E33.lat[i]]
            lon_sav2_ls90pm +=[DDS_E33.lon[j]]
        
for i in range(0,len(DDS_E33.lat)):
    for j in range(0,len(DDS_E33.lon)):
        tmp_ls270 = tmp_SY1[ls270,i,j]
        if np.any(tmp_ls270 < 2):
            continue
        else:
            lat_sav2_ls270pm+= [DDS_E33.lat[i]]
            lon_sav2_ls270pm +=[DDS_E33.lon[j]]

# SEASONAL DIURNAL AVERAGE
fig, [(ax1,ax2),(ax3,ax4)] = py.subplots(2,2,sharey=True, sharex=True,figsize=(15,8))
py.subplots_adjust(wspace=0.05)
fig1 = ax1.contourf(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf, 20,cmap=py.cm.viridis)
ax1.set_title('Proposed Landing Sites, Ls=90')
ax1.set_xlabel('Longitude')
ax1.xaxis.set_major_locator(MultipleLocator(60))
ax1.yaxis.set_major_locator(MultipleLocator(30))

# add landing sites
ax1.plot(longitude, latitude, 'ok',label='ROIs')
ax1.plot(lon_sav2_ls90, lat_sav2_ls90, 'ok', markersize=10,fillstyle='none', label='>2.2kW')
ax1.plot(lon_sav15_ls90, lat_sav15_ls90, 'vk', markersize=15, label='>14.2kW')
ax1.plot(lon_sav24_ls90, lat_sav24_ls90, 'sk', markersize=15, markerfacecolor='red', label='>24kW')
ax1.legend(loc='lower right')


# LS=270
fig2 = ax2.contourf(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf, 20,cmap=py.cm.viridis)
ax2.set_title('Proposed Landing Sites, Ls=270')


# add landing sites
ax2.plot(longitude, latitude, 'ok',label='ROIs')
ax2.plot(lon_sav2_ls270, lat_sav2_ls270, 'ok', markersize=10,fillstyle='none', label='>2.2kW')
ax2.plot(lon_sav15_ls270, lat_sav15_ls270, 'vk', markersize=15, label='>14.2kW')
ax2.plot(lon_sav24_ls270, lat_sav24_ls270, 'sk', markersize=15, markerfacecolor='red', label='>24kW')
ax2.legend(loc='lower right')

# NIGHTTIME SEASONAL AVERAGE

fig3 = ax3.contourf(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf, 20,cmap=py.cm.viridis)
ax3.set_title('Proposed Landing Sites, Ls=90')
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')

# add landing sites
ax3.plot(longitude, latitude, 'ok',label='ROIs')
ax3.plot(lon_sav2_ls90pm, lat_sav2_ls90pm, 'ok', markersize=10,fillstyle='none', label='>2.2kW')
ax3.plot(lon_sav15_ls90pm, lat_sav15_ls90pm, 'vk', markersize=15, label='>14.2kW')
ax3.plot(lon_sav24_ls90pm, lat_sav24_ls90pm, 'sk', markersize=15, markerfacecolor='red', label='>24kW')
ax3.legend(loc='lower right')


# LS=270
fig4 = ax4.contourf(DDS_SY1.lon, DDS_SY1.lat, DS_fixed.zsurf, 20,cmap=py.cm.viridis)
ax4.set_title('Proposed Landing Sites, Ls=270')
ax4.set_xlabel('Longitude')

# add landing sites
ax4.plot(longitude, latitude, 'ok',label='ROIs')
ax4.plot(lon_sav2_ls270pm, lat_sav2_ls270pm, 'ok', markersize=10,fillstyle='none', label='>2.2kW')
ax4.plot(lon_sav15_ls270pm, lat_sav15_ls270pm, 'vk', markersize=15, label='>14.2kW')
ax4.plot(lon_sav24_ls270pm, lat_sav24_ls270pm, 'sk', markersize=15, markerfacecolor='red', label='>24kW')
ax4.legend(loc='lower right')

py.savefig(f'{PATH}/WindEnergy_HIGHRES_FigS{ct}.eps',dpi=200)
