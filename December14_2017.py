#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:30:51 2020

@author: noahbrauer
"""




import h5py 
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles
from pyproj import Proj

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
import pyart


file = '2A.GPM.DPR.V8-20180723.20171214-S024517-E041750.021554.V06A.HDF5'


#Read in file, extract lat, lons, PSD

DPR = h5py.File(file, 'r')

#nscan is the nummber of scans in each granule 
#nray is the number of angle bins in each  scan; think of these as footprint scans (5 km in diameter for each footprint)
#nbin is the number of range bins in each ray (angle bins)
#nDSD: Parameters are N0 (number concentration) and D0 (mean drop diameter)


lat = DPR['NS']['Latitude'][:,:]    
lon = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:] #nscan x nray (7934,49)
za_ka = DPR['HS']['SLV']['zFactorCorrectedNearSurface'][:]
dsd = DPR['NS']['SLV']['paramDSD'][:] #nscan x nray x nbin x DSDmoment (7934,49,176,2)
freezing = DPR['NS']['VER']['heightZeroDeg'][:]
lwp = DPR['NS']['SLV']['precipWaterIntegrated'][:,:,0] #liquid water path
iwp = DPR['NS']['SLV']['precipWaterIntegrated'][:,:,1] #ice water path
precip_type = DPR['NS']['CSF']['typePrecip'][:]/10**7 #Convective vs stratiform 
precip_rate = precip = DPR['NS']['SLV']['precipRateNearSurface'][:]

#Extract Ka-band reflectivity

ka = DPR['HS']['SLV']['zFactorCorrected'][:]


#Add NaNs for zero values of iwp and lwp

lwp[lwp==0] = np.nan
iwp[iwp==0] = np.nan
precip_rate[precip_rate<=0] = np.nan



lat_ka = DPR['HS']['Latitude'][:]
lon_ka = DPR['HS']['Longitude'][:]

snow_rate = np.sqrt(z/120.)

snow_rate_ka = np.sqrt(za_ka/120)
snow_rate_ka[snow_rate_ka<=0] = np.nan

snow_rate[snow_rate<=0] = np.nan


#Compute the ratio of liquid water path to ice water path

ratio = lwp/iwp





###Define cross-section length

ind1 = np.where((lon[:,0]>=-83.75)) #Where are the DPR lons >= -100
ind2 = np.where((lon[:,0])<=-79.5) #Where are lons <= -85
ind3 = np.intersect1d(ind1,ind2) #Both conditions need to be satisfied here


#Change the ray to change cross-section location (i.e. the "27" in this case)
#%%

###Setup to 2D grid for plotting

x = 2.* 17 #48 degrees (from -17 to 17)
re = 6378. #radius of the earth
theta = -1*(x/2.) + (x/48.)*np.arange(0,49) #Split into equal degrees (from -17 to 17)
theta2  = np.ones(theta.shape[0]+1)*np.nan #Define an empty array (NaNs) with same shape as ray dimension
theta = theta - 0.70833333/2. #Shift degrees for plotting pruposes
theta2[:-1] = theta #remove last array dimension in theta so python won't get angry about shape
theta2[-1] = theta[-1] + 0.70833333
theta = theta2*(np.pi/180.) #Convert from degrees to radians

prh = np.ones((177,50))*np.nan #Define empty grid

for i in range(prh.shape[0]): #Loop over each range gate
    for j in range(prh.shape[1]): #Loop over each scan
            a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #Orbit height of 407 km 
            
            prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a)


prh_ka = np.ones((175,24))*np.nan
for i in range(prh_ka.shape[0]): #Loop over each range gate
    for j in range(prh_ka.shape[1]): #Loop over each scan
            a_ka = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #Orbit height of 407 km 
            
            prh_ka[i,j] = (176-(i))*0.125*np.cos(theta[j]+a_ka)




            
h2 = prh #where prh is the (range bins,ray)-space
h3 =np.ones((177,50))*np.nan

for i in range(h3.shape[1]):
    h3[:,i] = h2[::-1,i] #This reverses the vertical dimension so as indices increase, height increases


#%%
    
ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:] #Read in ku-band reflectivity; nscan x nray (554,49,176)
n0 = dsd[ind3,:,:,0]/10 #Read in the number concentration
d0 = dsd[ind3,:,:,1] #Read in the mean drop diameter  #Both have dimensions nscan x nray x nbin (554,49,176)
zeroDeg = freezing[ind3,:]


#Choose point along cross-section. Compute cross-track cross-sections

point = 44

ku_xtrack = ku[point,:,:]
n0_xtrack = n0[point,:,:]
d0_xtrack = d0[point,:,:]





#Cut all parameters so they are at same ray as above


ray = 20

ku = ku[:,ray,:]
n0 = n0[:,ray,:]
d0 = d0[:,ray,:]
zero_deg_isotherm = zeroDeg[:,ray]/1000 #Convert from meters to km

#Take lats and lons along same ray
lons = DPR['NS']['Longitude'][ind3,ray]
lats = DPR['NS']['Latitude'][ind3,ray]

#Compute cross-track cross-sections








#Choose a starting point, then calculate distance
lat0 = lats[0]
lon0 = lons[0]


p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0) #Define a projection and set starting lat an lon to same point as above

#Define a 2D array for plotting purposes

lat_3d = np.ones(ku.shape)*np.nan
lon_3d = np.ones(ku.shape)*np.nan

for i in range(ku.shape[0]):
    lat_3d[i,:] = lats[i] 
    lon_3d[i,:] = lons[i]  
        

x,y = p(lon_3d,lat_3d) #Now convert degrees to distance (in km)
R_gpm = np.sqrt(x**2 + y**2)*np.sign(x) #Keeps sign of number; converts to radial distance 







#%%

#Now do this for the cross-track
lons_xtrack = DPR['NS']['Longitude'][point,:]
lats_xtrack = DPR['NS']['Latitude'][point,:]

#Compute cross-track cross-sections


#Choose a starting point, then calculate distance
lat0_xtrack = lats_xtrack[0]
lon0_xtrack = lons_xtrack[0]


p_xtrack = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0) #Define a projection and set starting lat an lon to same point as above

#Define a 2D array for plotting purposes

lat_3d_xtrack = np.ones(ku_xtrack.shape)*np.nan
lon_3d_xtrack = np.ones(ku_xtrack.shape)*np.nan

for i in range(ku_xtrack.shape[0]):
    lat_3d_xtrack[i,:] = lats_xtrack[i] 
    lon_3d_xtrack[i,:] = lons_xtrack[i]  
        




y_xtrack = np.ones((ku_xtrack.shape))

for j in range(ku_xtrack.shape[0]):
    
    lat_3d_xtrack[i,:] = lats_xtrack[i]
    lon_3d_xtrack[i,:] = lons_xtrack[i]
    
    
x_xtrack,y_xtrack = p_xtrack(lon_3d_xtrack,lat_3d_xtrack) #Now convert degrees to distance (in km)
R_gpm_xtrack = np.sqrt(x_xtrack**2 + y_xtrack**2)*np.sign(x_xtrack) #Keeps sign of number; converts to radial distance    

#%%


#Reverse range gate order for all parameters

ku = ku[:,::-1]
n0 = n0[:,::-1]
d0 = d0[:,::-1]


ku_xtrack = ku_xtrack[:,::-1]
n0_xtrack = n0_xtrack[:,::-1]
d0_xtrack = d0_xtrack[:,::-1]



ku = np.ma.masked_where(ku<=12, ku)
ku_xtrack = np.ma.masked_where(ku_xtrack<=12, ku_xtrack) #Mask all the bad points in ku data
y = np.ones([ku.shape[0], ku.shape[1]]) #Define an empty array

#Define the appropriate range bins
h4 = h3[:,ray] #This number MUST match the same ray being used
for i in range(y.shape[1]):
    y[:,i] = h4[i]
    



    
    
##Now compute 2D histogram for ku-band reflectivity moving along track
#%%    

   
def create_histogram(input_values, num_bins, min_value, max_value):
    """Creates a histogram with uniform bin-spacing.
    N = number of input values
    K = number of bins
    :param input_values: length-N numpy array of input values (to be binned).
    :param num_bins: Number of bins.
    :param min_value: Minimum value to include in histogram.  Any input value <
        `min_value` will be assigned to the first bin.
    :param max_value: Maximum value to include in histogram.  Any input value >
        `max_value` will be assigned to the last bin.
    :return: input_to_bin_indices: length-N numpy array of bin indices.  If
        input_values[i] = j, the [i]th input value belongs in the [j]th bin.
    :return: num_examples_by_bin: length-K numpy array, where the [j]th value is
        the number of inputs assigned to the [j]th bin.
    """

    bin_cutoffs = np.linspace(min_value, max_value, num=num_bins + 1)
    input_to_bin_indices = np.digitize(
        input_values, bin_cutoffs, right=False) - 1
    input_to_bin_indices[input_to_bin_indices < 0] = 0
    input_to_bin_indices[input_to_bin_indices > num_bins - 1] = num_bins - 1

    num_examples_by_bin = np.full(num_bins, -1, dtype=int)
    for j in range(num_bins):
        num_examples_by_bin[j] = np.sum(input_to_bin_indices == j)

    return input_to_bin_indices, num_examples_by_bin     


ku_hist = np.zeros((40,176))


for i in range(ku.shape[1]):
    
    ku_hist[:,i] = create_histogram(ku[:,i], 40,10,50)[1]    
    
#Take the transpose so it's row = bins, column = height
    
ku_hist_no_nan = ku_hist.copy() 
 

#NaN out zeros:

ku_hist[ku_hist==0] = np.nan
ku_hist[ku_hist>30] = np.nan


x_bins = np.arange(10,50,1)
y_bins = np.arange(0,21.988032,0.124932)

fig,ax = plt.subplots(figsize=(14,14)) 
plt.pcolormesh(x_bins,y_bins,ku_hist.T, cmap = 'pyart_NWSRef')

xtick_label_size = 17
ytick_label_size = 17
tick_label_size = 24
title_size = 28


#Plot the mean height of the zero degree isotherm


zero_deg_mean = np.nanmean(zero_deg_isotherm)

#pm2 = plt.plot(y_bins,zero_deg_mean, '--', color = 'k')
pm2 = plt.hlines(zero_deg_mean,x_bins[0],x_bins[9], linestyles = 'dashed', color = 'k')


xlabels = np.arange(12.5,65,5)
ylabels = np.arange(0,13.5,0.5)

plt.xticks(xlabels)
plt.yticks(ylabels)


plt.xlim(10,40)
plt.ylim(0,8)
plt.clim(0,30)

plt.ylabel('Altitude (km)', size = 16)
plt.xlabel('Bin range (dBZ)', size = 16)

ax.set_xticklabels(['<10','[10,15)','[15,20)', '[20,25)','[25,30)','[30,35)','[35,40)','[40,45)','[45,50)','[50,55)'], size = xtick_label_size)
plt.yticks(size =ytick_label_size )        
plt.title('Cleveland, OH 12/14/2017 Along-Track KuPR', size = title_size)


cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
cbar.set_label(label = '[Frequency]',size = tick_label_size)

plt.show()   









#%%%
'''

print(np.nanmin(ku))
print(np.nanmax(ku))
#Okay, now lets try and plot as a line plot
#Ku data is organized as along-track distance by height

y_bins_vert = np.arange(0,43.75,step = 0.25)

ku_zeros = np.array(ku, dtype = np.float)
ku_zeros[ku_zeros<=0] = np.nan

ku_25 = np.nanpercentile(ku_zeros,25,axis = 1)
ku_50 = np.nanpercentile(ku_zeros,50,axis = 1)
ku_75 = np.nanpercentile(ku_zeros,75, axis = 1)
ku_95 = np.nanpercentile(ku_zeros,95, axis = 1)


fig,ax = plt.subplots(figsize=(14,14)) 
plt.plot(ku_25,y_bins_vert, label = '25q')
plt.plot(ku_50,y_bins_vert, label = '50q')
plt.plot(ku_75,y_bins_vert, label = '75q')
plt.plot(ku_95,y_bins_vert, label = '95q')





xlabels = np.arange(12.5,40,2.5)
ylabels = np.arange(0,13.5,0.5)

plt.xticks(xlabels)
plt.yticks(ylabels, size = xtick_label_size)



plt.ylim(0,13)


plt.ylabel('Altitude (km)', size = 16)
plt.xlabel('dBZ', size = 16)

ax.set_xticklabels(['12.5','15','17.5','20','22.5','25','27.5','30','32.5','35','37.5','40'], size = xtick_label_size)

        
plt.title('Buffalo, NY 3/2 Along-Track KuPR', size = title_size)
plt.legend()
plt.show()
'''







        



###Calculate echo top height:


#echo_top_height_index = np.array([15,34,34,27,16,17,17,20,28,34,38,39,44,40,44,34,40,36,47,37,29,34,26,36,37,41,39,42,41,40,37,40,37,39,40,36,37,36,34,33,35,37,34,33,32,33,35,33,33,35,37,28,27,31,32,35,31,31,28,30,40,39,43,40,42,37,32,27,26,33,30,28,31,29,34,28,29,26,23,27,])

echo_top_height_index = np.array([16,18,17,15,18,17,15,19,20,15,19,25,22,25,23,21,18,18,22,22,20,20,14,14,18,15,16,17,14,13,15,17,18,])
echo_top_height = echo_top_height_index*0.124932

echo_top_mean = np.nanmean(echo_top_height)
echo_top_std = np.nanstd(echo_top_height)

bins = np.arange(1,8,step = 0.5)
plt.hist(echo_top_height,bins = bins,edgecolor='black', linewidth=1.2, color = 'blue')
plt.xlabel('Echo Top Height (km)', fontsize = 20)
plt.ylabel('Frequency', fontsize = 20)
plt.title('12/14/2017 Along-Track Echo Top Height', fontsize = 18)
frequency_bins = np.arange(0,35, step = 5)
plt.axvline(echo_top_mean, color='k', linestyle='dashed', linewidth=4)
plt.xticks(size = 16)
plt.yticks(frequency_bins,size = 16)


print(np.nanmin(echo_top_height))




 
    





    







#%%    


#Remove the values less than or equal to zero

n0[n0<=0] = np.nan
d0[d0<=0] = np.nan

n0_xtrack[n0_xtrack<=0] = np.nan
d0_xtrack[d0_xtrack<=0] = np.nan

    
    
#%% 
    
#Now we plot! #N-S (along track first)  


#Plot mean drop size first

dist = np.array([0,300])
altitude = np.array([0,10])

x_loc = 500
y_loc = 14.5
label = '@NOAABrauer'

plt.figure(figsize=(10,10))

vmax = 3
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()

label_size = 20

cmin =0.; cmax = 3.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, d0, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 12/14/2017 $D_{M}$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(altitude[0],altitude[1])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[mm]',size = label_size)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.clim(0,3)


#plt.text(x_loc, y_loc, label)


plt.show()    
    


####Number concentration (liquid water content)

plt.figure(figsize=(10,10))

vmax = 6
vmin = 1

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =1.; cmax = 6.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, n0, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 12/14/2017 $log_{10 }(N_{W})$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(altitude[0],altitude[1])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[mm $m^{-3}$]',size = label_size)
plt.clim(1,6)
#plt.text(x_loc, y_loc, label)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.show()    
    

###And lastly Ku-band


plt.figure(figsize=(10,10))

vmax = 60
vmin =12

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 12/14/2017 KuPR ', size = 20)
#plt.xlim(300,450)
plt.xlim(dist[0], dist[1])
plt.ylim(altitude[0],altitude[1])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.clim(12,60)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

#plt.text(x_loc, y_loc, label)
plt.show()   


#%%
'''
#Now for cross-track

dist = np.array([0,500])
altitude = np.array([0,10])

x_loc = 500
y_loc = 14.5
label = '@NOAABrauer'

plt.figure(figsize=(10,10))

vmax = 3
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()

label_size = 20

cmin =0.; cmax = 3.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm_xtrack/1000., y, d0_xtrack.T, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm_xtrack/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 3/8 1452 UTC $D_{M}$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(altitude[0],altitude[1])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[mm]',size = label_size)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.clim(0,3)


#plt.text(x_loc, y_loc, label)


plt.show()    
    


####Number concentration (liquid water content)

plt.figure(figsize=(10,10))

vmax = 6
vmin = 1

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =1.; cmax = 6.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm_xtrack/1000., y, n0_xtrack, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm_xtrack/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 3/8 1452 UTC $log_{10 }(N_{W})$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(altitude[0],altitude[1])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[mm $m^{-3}$]',size = label_size)
plt.clim(1,6)
#plt.text(x_loc, y_loc, label)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.show()    
    

###And lastly Ku-band


plt.figure(figsize=(10,10))

vmax = 60
vmin =12

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm_xtrack/1000., y, ku_xtrack, cmap='pyart_NWSRef',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm_xtrack/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 3/8 1452 UTC KuPR ', size = 20)
#plt.xlim(300,450)
plt.xlim(dist[0], dist[1])
plt.ylim(altitude[0],altitude[1])
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.clim(12,60)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

#plt.text(x_loc, y_loc, label)
plt.show()  

'''



#%%

#Determine median value of ind3 (intersection between longitudes)
cross_track_index = int(len(ind3)/2)+ind3[0]



#Let's plot a map with the center point location for our cross-sections

#Mask near surface reflectivity values

z = np.ma.masked_where(z<=12, z)



cmin = 12.; cmax = 70.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-79]); ylim = np.array([39,43])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,z,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()

parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-80,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])

cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[dBZ]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)
   
plt.title('GPM Overpass 12/14/2017 KuPR', size = 20)




print(np.nanmin(ratio))
print(np.nanmax(ratio))
#%%


print(np.nanmax(precip_rate))

cmin = 20.; cmax = 50.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-79]); ylim = np.array([39,43])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,z,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-80,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[mm/hour]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)
   
plt.title('GPM Overpass 12/14/2017 Near Surface Precipitation Rate', size = 20)



#%%
#Now plot


#Set categories for precip type:

precip_type[precip_type<1] = np.nan



convective = precip_type.copy()
convective[convective<=1] = np.nan

print(np.nanmax(convective))
print(np.nanmin(convective))


text_size = 20



cmin = 0.; cmax = 150.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-79]); ylim = np.array([39,43])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,z,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)
m.contour(lon,lat,convective, color = 'black')


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-80,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[LWP/IWP]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)

 
plt.title('GPM Overpass 12/14/2017 LWP/IWP Ratio', size = 20)

print(np.nanmax(iwp))

#%%


cmin = 0.; cmax = 100.; cint = 5.; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-79]); ylim = np.array([39,43])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,iwp/10,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)



##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,ray],lat[cross_track_index,ray],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)


m.drawcounties()
parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-80,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[$gm^{-2}$]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)

 
plt.title('GPM Overpass 12/14/2017 IWP', size = 20)



#%%

colormap = ['white','green', 'yellow', 'red']
 
from matplotlib.colors import ListedColormap
cmap_precip = ListedColormap(colormap)     

print(np.nanmax(precip_type))



plt.figure(figsize=(10,10))

cmin = 1.; cmax = 4.; cint = 1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_precip,lut=nlevs)
xlim = np.array([-83,-79]); ylim = np.array([39,43])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,precip_type,clevs,cmap=cmap_precip,extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 20, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])

parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-80,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])

m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.drawcounties()
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_yticklabels(['None','Stratiform', 'Convective', 'Other'], size = 20)



plt.title('GPM Overpass 12/14/2017 Precipitation Category', size = 20)
#%%



cmin = 0.2; cmax = 0.9; cint = 0.025; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-79]); ylim = np.array([39,43])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,snow_rate,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)



##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,ray],lat[cross_track_index,ray],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)


m.drawcounties()
parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-84,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[mm/hour]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)

 
plt.title('GPM Overpass 12/14/2017 Snowfall Rate (KuPR)', size = 20)




#Now for KaPR


cmin = 0.2; cmax = 0.9; cint = 0.025; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-83,-79]); ylim = np.array([39,43])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon_ka,lat_ka,snow_rate_ka,clevs,cmap='pyart_NWSRef',extend='both')
cs2 = m.plot(lon_ka[:,0]+0.03,lat_ka[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon_ka[:,-1]-0.03,lat_ka[:,-1]+0.03,'--k',zorder=4)
#m.plot(lon_ka[ind3,ray],lat_ka[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


'''
##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,ray],lat[cross_track_index,ray],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
'''

m.drawcounties()
parallels = np.arange(30,50,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-84,-66, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[mm/hour]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)

 
plt.title('GPM Overpass 12/14/2017 Snowfall Rate (KaPR)', size = 20)

'''
latitude_path = [41.8,42.2,42.8,42.8]
longitude_path = [-72.3,-72.3,-71.8,-71.4]



x,y = m(longitude_path,latitude_path)
m.plot(x,y,'k-', linewidth = 4)
'''

