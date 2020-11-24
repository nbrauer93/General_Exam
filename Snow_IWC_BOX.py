#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 12:34:22 2020

@author: noahbrauer
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 10:20:23 2020

@author: noahbrauer
"""


import netCDF4
import gzip
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import pyart



import os
import shutil
import tempfile



def open_netcdf(fname):
    if fname.endswith(".gz"):
        infile = gzip.open(fname, 'rb')
        tmp = tempfile.NamedTemporaryFile(delete=False)
        shutil.copyfileobj(infile, tmp)
        infile.close()
        tmp.close()
        data = netCDF4.Dataset(tmp.name)
        os.unlink(tmp.name)
    else:
        data = netCDF4.Dataset(fname)
    return data

#Use Ka-band to measure with high sensitivity

file = 'GRtoDPR.KBOX.180308.22860.V06A.KU.NS.1_21.nc.gz'





nc = open_netcdf(file)



#For the Ku-band file
latitude = nc.variables['latitude'][:].T
longitude = nc.variables['longitude'][:].T
z = nc.variables['ZFactorCorrected'][:] #GPM DPR attenuation-corrected reflectivity
nw = nc.variables['Nw'][:]/10 #converts dBNw t0 long10(nw) from GPM DPR
nw_ground = nc.variables['GR_Nw'][:]/10
dm = nc.variables['Dm'][:] #From GPM DPR
dm_ground = nc.variables['GR_Dm'][:] #Ground-radar median volume diameter
rain_rate_ground = nc.variables['GR_RR_rainrate'][:]
rain_rate_DPR = nc.variables['PrecipRateSurface'][:] #DPR near surface precip rate
rain_rate_GMI = nc.variables['SurfPrecipTotRate'][:] #GMI rainfall rate
elevation_angle = nc.variables['elevationAngle'][:]
site_elev = nc.variables['site_elev'][:]*1000 #in meters
gr_z = nc.variables['GR_Z'][:]
zdr = nc.variables['GR_Zdr'][:]
kdp = nc.variables['GR_Kdp'][:]




#Now calculate S(Kdp,Z) and IWC(Kdp,Z):


s_kdp_z = 1.48*(kdp**(0.615))*(gr_z**(0.33)) #mm/hour
iwc_kdp_z = 0.71*(kdp**(0.66))*(gr_z**(0.28)) #gm^-3

print(np.nanmax(s_kdp_z))
print(np.nanmax(iwc_kdp_z))
print(np.nanmin(iwc_kdp_z))

#Extract lat-lon of DPR footprints


latitude_fp = nc.variables['DPRlatitude'][:]
longitude_fp = nc.variables['DPRlongitude'][:]


#Remove bad values

nw[nw<=0] = np.nan
dm[dm<=0] = np.nan
nw_ground[nw_ground<=0] = np.nan
dm_ground[dm_ground<=0] = np.nan
rain_rate_ground[rain_rate_ground<=0] = np.nan
rain_rate_DPR[rain_rate_DPR<=0] = np.nan
rain_rate_GMI[rain_rate_GMI<=0] = np.nan
gr_z[gr_z<=0] = np.nan
zdr[zdr<0] = np.nan
kdp[kdp<=0] = np.nan
s_kdp_z[s_kdp_z<=0] = np.nan
iwc_kdp_z[iwc_kdp_z<=0] = np.nan
z[z<=0] = np.nan






#Select the latitude and longitude of Houston: 




elev_angle = 0
z_lowest = z[elev_angle,:]
nw_lowest = nw[elev_angle,:]
nw_ground_lowest = nw_ground[elev_angle,:]
dm_lowest = dm[elev_angle,:]
dm_ground_lowest = dm_ground[elev_angle,:]
rain_rate_lowest_ground = rain_rate_ground[elev_angle,:]
s_kdp_z_lowest = s_kdp_z[elev_angle,:]
iwc_kdp_z_lowest = iwc_kdp_z[elev_angle,:]



snow_z = np.sqrt(gr_z[elev_angle,:]/120)
snow_ku = np.sqrt(z_lowest/120)



#Compute the dual-frequency ratio


#dfr = 10*np.log(ku_ka_lowest) - 10*np.log(z_ka_lowest)



#Now let's compare quantities

#Compare the fraction of estimated rainfall rate from ground radar to near-surface rainfall rate of GPM DPR:

rainfall_rate_fracion = rain_rate_lowest_ground/rain_rate_DPR



#Determine beam height and distance from the radar; Equations from Doviak and Zrnic

deg_to_rad = np.pi/180. #Convert from degrees to radians
a_e = (8494.66667)*1000
k = 4/3 #Use small angle approximation

#Compute distance from the radar location

radar_lat = nc.variables['site_lat'][:]
radar_lon = nc.variables['site_lon'][:]

#Now moving away from the radar site, calculate the difference in degrees in both zonal and meridional directions:




distance_x_deg = np.ones(longitude.shape)*np.nan
distance_y_deg = np.ones(latitude.shape)*np.nan

for i in range(latitude.shape[0]):
    for j in range(latitude.shape[1]):
        
        distance_x_deg[i,j] = radar_lon - longitude[i,j]
        distance_y_deg[i,j] = radar_lat - latitude[i,j]
        
        

#Convert grid to utm coordinates

pp = Proj(proj = 'utm', zone = 10, ellps = 'WGS84', preserve_units = False)


        
'''         
#Now compute degrees to distance from the radar:
        
distance_x = distance_x_deg*111000 #In meters
distance_y = distance_y_deg*111000         
'''
#And compute the total distance of each point (s):

#s = np.sqrt((distance_x**2) + (distance_y**2))

#Compute the sine (in radians) of each elevation angle:

elev_angle_rad = np.sin(elevation_angle*deg_to_rad)



#Now lets compute beam height for each elevation angle using equation from Doviak and Zrnic:
'''
beam_height = np.ones(s.shape)*np.nan

for k in range(0,len(elevation_angle)):
    
    beam_height = site_elev - k*a_e + np.sqrt((s**2) + ((k**2)*(a_e**2) + 2*s*k*a_e*np.sin(elev_angle_rad[k])))


#Convert beam height from m to km
    
beam_height_m = beam_height/1000    

plt.figure(figsize=(18,10))
plt.imshow(z[::-1,:])
    
'''  
#Compute the difference between Zh and KuPR to determine bias:


z_bias = gr_z[elev_angle,:] - z_lowest      

        
    
        



    
    
    






#%%

#from matplotlib.colors import ListedColormap

colormap = ['white','dodgerblue', 'deepskyblue', 'lawngreen', 'lightgreen', 'green', 'gold', 'darkorange', 'red', 'firebrick']

z_color = np.empty(667,dtype = 'str')
z_color = []

#for i in range(len(z_color)):
for i in range(len(z_lowest)):

    if z_lowest[i]<=5:
        #z_color[i] = colormap[0]
        z_color.append(colormap[0])
    
    elif z_lowest[i]>5 and z_lowest[i]<=10:
        #z_color[i] = colormap[1]
        z_color.append(colormap[1])
    
    elif z_lowest[i]>10 and z_lowest[i]<=15:
        #z_color[i] = colormap[2]
        z_color.append(colormap[2])
    
    elif z_lowest[i]>15 and z_lowest[i]<=20:
        #z_color[i] = colormap[3]
        z_color.append(colormap[3])
    
    elif z_lowest[i]>20 and z_lowest[i]<=25:
        #z_color[i] = colormap[4]
        z_color.append(colormap[4])
    
    elif z_lowest[i]>25 and z_lowest[i]<=30:
        #z_color[i] = colormap[5]
        z_color.append(colormap[5])
        
    elif z_lowest[i]>30 and z_lowest[i]<=35:
        #z_color[i] = colormap[6]
        z_color.append(colormap[6])
    
    elif z_lowest[i]>35 and z_lowest[i]<=40:
        #z_color[i] = colormap[7]
        z_color.append(colormap[7])
        
    elif z_lowest[i]>40 and z_lowest[i]<=45:
        #z_color[i] = colormap[8]
        z_color.append(colormap[8])
        
    elif z_lowest[i] == -100:
        #z_color[i] = colormap[0]
        z_color.append(colormap[0])
        
        
from matplotlib.colors import ListedColormap
cmap_z = ListedColormap(colormap)        
    

print(np.nanmax(nw))
print(np.nanmin(nw))
#Now do the same thing for log10(Nw) 


nw_color = []

for i in range(len(nw)):

    if nw_lowest[i]<=1:
        
        nw_color.append(colormap[0])
    
    elif nw_lowest[i]>1 and nw_lowest[i]<=1.5:
        
        nw_color.append(colormap[1])
    
    elif nw_lowest[i]>1.5 and nw_lowest[i]<=2:
        
        nw_color.append(colormap[2])
    
    elif nw_lowest[i]>2 and nw_lowest[i]<=2.5:
        
        nw_color.append(colormap[3])
    
    elif nw_lowest[i]>2.5 and nw_lowest[i]<=3:
        
        nw_color.append(colormap[4])
    
    elif nw_lowest[i]>3 and nw_lowest[i]<=3.5:
        
        nw_color.append(colormap[5])
        
    elif nw_lowest[i]>3.5 and nw_lowest[i]<=4:
        
        nw_color.append(colormap[6])
    
    elif nw_lowest[i]>4 and nw_lowest[i]<=4.5:
       
        nw_color.append(colormap[7])
        
    elif nw_lowest[i]>4.5 and nw_lowest[i]<=5:
        
        nw_color.append(colormap[8])
        
    
#Now for dm
        
        
print(np.nanmin(dm))
print(np.nanmax(dm))  

dm_colormap = ['white','dodgerblue', 'deepskyblue', 'lawngreen', 'lightgreen', 'green', 'yellow', 'gold', 'darkorange', 'red', 'firebrick', 'blueviolet', 'magenta', 'deeppink', 'hotpink', 'pink']      
        
dm_color = []    

for i in range(len(dm)):

    if dm_lowest[i]<=0.5:
        
        dm_color.append(dm_colormap[0])
    
    elif dm_lowest[i]>0.5 and dm_lowest[i]<=0.75:
        
        dm_color.append(dm_colormap[1])
    
    elif dm_lowest[i]>0.75 and dm_lowest[i]<=1.0:
        
        dm_color.append(dm_colormap[2])
    
    elif dm_lowest[i]>1.0 and dm_lowest[i]<=1.25:
        
        dm_color.append(dm_colormap[3])
    
    elif dm_lowest[i]>1.25 and dm_lowest[i]<=1.5:
        
        dm_color.append(dm_colormap[4])
    
    elif dm_lowest[i]>1.5 and dm_lowest[i]<=1.75:
        
        dm_color.append(dm_colormap[5])
        
    elif dm_lowest[i]>1.75 and dm_lowest[i]<=2.0:
        
        dm_color.append(dm_colormap[6])
    
    elif dm_lowest[i]>2.0 and dm_lowest[i]<=2.25:
       
        dm_color.append(dm_colormap[7])
        
    elif dm_lowest[i]>2.25 and dm_lowest[i]<=2.5:
        
        dm_color.append(dm_colormap[8])
        
    elif dm_lowest[i]>2.5 and dm_lowest[i]<=2.75:
        
        dm_color.append(dm_colormap[9])
        
    elif dm_lowest[i]>2.75 and dm_lowest[i]<=3.0:
        
        dm_color.append(dm_colormap[10])
        
    elif dm_lowest[i]>3.0 and dm_lowest[i]<=3.25:
        
        dm_color.append(dm_colormap[11])
        
    elif dm_lowest[i]>3.25 and dm_lowest[i]<=3.5:
        
        dm_color.append(dm_colormap[12])
        
    elif dm_lowest[i]>3.5 and dm_lowest[i]<=3.75:
        
        dm_color.append(dm_colormap[13]) 
        
    elif dm_lowest[i]>3.75 and dm_lowest[i]<=4.0:
        
        dm_color.append(dm_colormap[14])
            
            
cmap_dm = ListedColormap(dm_colormap)           


#Now for rainfall rate:


rainfall_color = []

#for i in range(len(z_color)):
for i in range(len(rain_rate_DPR)):

    if rain_rate_DPR[i]<=1:
        #z_color[i] = colormap[0]
        rainfall_color.append(dm_colormap[0])
    
    elif rain_rate_DPR[i]>1 and rain_rate_DPR[i]<=5:
        #z_color[i] = colormap[1]
        rainfall_color.append(dm_colormap[1])
    
    elif rain_rate_DPR[i]>5 and rain_rate_DPR[i]<=10:
        #z_color[i] = colormap[2]
        rainfall_color.append(dm_colormap[2])
    
    elif rain_rate_DPR[i]>10 and rain_rate_DPR[i]<=15:
        #z_color[i] = colormap[3]
        rainfall_color.append(dm_colormap[3])
    
    elif rain_rate_DPR[i]>15 and rain_rate_DPR[i]<=20:
        #z_color[i] = colormap[4]
        rainfall_color.append(dm_colormap[4])
    
    elif rain_rate_DPR[i]>20 and rain_rate_DPR[i]<=25:
        #z_color[i] = colormap[5]
        rainfall_color.append(dm_colormap[5])
        
    elif rain_rate_DPR[i]>25 and rain_rate_DPR[i]<=30:
        #z_color[i] = colormap[6]
        rainfall_color.append(dm_colormap[6])
    
    elif rain_rate_DPR[i]>30 and rain_rate_DPR[i]<=35:
        #z_color[i] = colormap[7]
        rainfall_color.append(dm_colormap[7])
        
    elif rain_rate_DPR[i]>35 and rain_rate_DPR[i]<=40:
        #z_color[i] = colormap[8]
        rainfall_color.append(dm_colormap[8])
        
    elif rain_rate_DPR[i]>40 and rain_rate_DPR[i]<=45:
        
        rainfall_color.append(dm_colormap[9])
        
    elif rain_rate_DPR[i]>45 and rain_rate_DPR[i]<=50:
        
        rainfall_color.append(dm_colormap[10])
        
    elif rain_rate_DPR[i]>50 and rain_rate_DPR[i]<=55:
        
        rainfall_color.append(dm_colormap[11])
        
    elif rain_rate_DPR[i]>55 and rain_rate_DPR[i]<=60:
        
        rainfall_color.append(dm_colormap[11])



#Snowfall rate:
        
        
cmap_snow_array = np.array([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.5,2.,2.5,3.,3.5,4.,4.5])  


snowfall_color = []



for i in range(len(s_kdp_z_lowest)):

    if s_kdp_z_lowest[i]<=0.3:
        #z_color[i] = colormap[0]
        snowfall_color.append(dm_colormap[0])
    
    elif s_kdp_z_lowest[i]>0.3 and s_kdp_z_lowest[i]<=0.4:
        #z_color[i] = colormap[1]
        snowfall_color.append(dm_colormap[1])
    
    elif s_kdp_z_lowest[i]>0.4 and s_kdp_z_lowest[i]<=0.5:
        #z_color[i] = colormap[2]
        snowfall_color.append(dm_colormap[2])
    
    elif s_kdp_z_lowest[i]>0.5 and s_kdp_z_lowest[i]<=0.6:
        #z_color[i] = colormap[3]
        snowfall_color.append(dm_colormap[3])
    
    elif s_kdp_z_lowest[i]>0.6 and s_kdp_z_lowest[i]<=0.7:
        #z_color[i] = colormap[4]
        snowfall_color.append(dm_colormap[4])
    
    elif s_kdp_z_lowest[i]>0.7 and s_kdp_z_lowest[i]<=0.8:
        #z_color[i] = colormap[5]
        snowfall_color.append(dm_colormap[5])
        
    elif s_kdp_z_lowest[i]>0.8 and s_kdp_z_lowest[i]<=0.9:
        #z_color[i] = colormap[6]
        snowfall_color.append(dm_colormap[6])
    
    elif s_kdp_z_lowest[i]>0.9 and s_kdp_z_lowest[i]<=1.:
        #z_color[i] = colormap[7]
        snowfall_color.append(dm_colormap[7])
        
    elif s_kdp_z_lowest[i]>1. and s_kdp_z_lowest[i]<=1.5:
        #z_color[i] = colormap[8]
        snowfall_color.append(dm_colormap[8])
        
    elif s_kdp_z_lowest[i]>1.5 and s_kdp_z_lowest[i]<=2.:
        
        snowfall_color.append(dm_colormap[9])
        
    elif s_kdp_z_lowest[i]>2. and s_kdp_z_lowest[i]<=2.5:
        
        snowfall_color.append(dm_colormap[10])
        
    elif s_kdp_z_lowest[i]>2.5 and s_kdp_z_lowest[i]<=3.:
        
        snowfall_color.append(dm_colormap[11])
        
    elif s_kdp_z_lowest[i]>3. and s_kdp_z_lowest[i]<=3.5:
        
        snowfall_color.append(dm_colormap[12])       
        
    elif s_kdp_z_lowest[i]>3.5 and s_kdp_z_lowest[i]<=4.:
        
        snowfall_color.append(dm_colormap[13])
        
    elif s_kdp_z_lowest[i]>4. and s_kdp_z_lowest[i]<=4.5:
        
        snowfall_color.append(dm_colormap[14])



        
iwc_color = []



for i in range(len(iwc_kdp_z_lowest)):

    if iwc_kdp_z_lowest[i]<=0.1:
        #z_color[i] = colormap[0]
        iwc_color.append(dm_colormap[0])
    
    elif iwc_kdp_z_lowest[i]>0.1 and iwc_kdp_z_lowest[i]<=0.2:
        #z_color[i] = colormap[1]
        iwc_color.append(dm_colormap[1])
    
    elif iwc_kdp_z_lowest[i]>0.2 and iwc_kdp_z_lowest[i]<=0.3:
        #z_color[i] = colormap[2]
        iwc_color.append(dm_colormap[2])
    
    elif iwc_kdp_z_lowest[i]>0.3 and iwc_kdp_z_lowest[i]<=0.4:
        #z_color[i] = colormap[3]
        iwc_color.append(dm_colormap[3])
    
    elif iwc_kdp_z_lowest[i]>0.4 and iwc_kdp_z_lowest[i]<=0.5:
        #z_color[i] = colormap[4]
        iwc_color.append(dm_colormap[4])
    
    elif iwc_kdp_z_lowest[i]>0.5 and iwc_kdp_z_lowest[i]<=0.6:
        #z_color[i] = colormap[5]
        iwc_color.append(dm_colormap[5])
        
    elif iwc_kdp_z_lowest[i]>0.6 and iwc_kdp_z_lowest[i]<=0.7:
        #z_color[i] = colormap[6]
        iwc_color.append(dm_colormap[6])
    
    elif iwc_kdp_z_lowest[i]>0.7 and iwc_kdp_z_lowest[i]<=0.8:
        #z_color[i] = colormap[7]
        iwc_color.append(dm_colormap[7])
        
    elif iwc_kdp_z_lowest[i]>0.8 and iwc_kdp_z_lowest[i]<=0.9:
        #z_color[i] = colormap[8]
        iwc_color.append(dm_colormap[8])
        
    elif iwc_kdp_z_lowest[i]>0.9 and iwc_kdp_z_lowest[i]<=1.0:
        
        iwc_color.append(dm_colormap[9])
        
    elif iwc_kdp_z_lowest[i]>1.0 and iwc_kdp_z_lowest[i]<=1.1:
        
        iwc_color.append(dm_colormap[10])
        
    elif iwc_kdp_z_lowest[i]>1.1 and iwc_kdp_z_lowest[i]<=1.2:
        
        iwc_color.append(dm_colormap[11])
        
    elif iwc_kdp_z_lowest[i]>1.2 and iwc_kdp_z_lowest[i]<=1.3:
        
        iwc_color.append(dm_colormap[12])       
        
    elif iwc_kdp_z_lowest[i]>1.3 and iwc_kdp_z_lowest[i]<=1.4:
        
        iwc_color.append(dm_colormap[13])
        
    elif iwc_kdp_z_lowest[i]>1.4 and iwc_kdp_z_lowest[i]<=1.5:
        
        iwc_color.append(dm_colormap[14])
        
#Now for dual-frequency ratio


#Now for S(Z):
        
snowfall_color_z = []



for i in range(len(snow_z)):

    if snow_z[i]<=0.3:
        #z_color[i] = colormap[0]
        snowfall_color_z.append(dm_colormap[0])
    
    elif snow_z[i]>0.3 and snow_z[i]<=0.4:
        #z_color[i] = colormap[1]
        snowfall_color_z.append(dm_colormap[1])
    
    elif snow_z[i]>0.4 and snow_z[i]<=0.5:
        #z_color[i] = colormap[2]
        snowfall_color_z.append(dm_colormap[2])
    
    elif snow_z[i]>0.5 and snow_z[i]<=0.6:
        #z_color[i] = colormap[3]
        snowfall_color_z.append(dm_colormap[3])
    
    elif snow_z[i]>0.6 and snow_z[i]<=0.7:
        #z_color[i] = colormap[4]
        snowfall_color_z.append(dm_colormap[4])
    
    elif snow_z[i]>0.7 and snow_z[i]<=0.8:
        #z_color[i] = colormap[5]
        snowfall_color_z.append(dm_colormap[5])
        
    elif snow_z[i]>0.8 and snow_z[i]<=0.9:
        #z_color[i] = colormap[6]
        snowfall_color_z.append(dm_colormap[6])
    
    elif snow_z[i]>0.9 and snow_z[i]<=1.:
        #z_color[i] = colormap[7]
        snowfall_color_z.append(dm_colormap[7])
        
    elif snow_z[i]>1. and snow_z[i]<=1.5:
        #z_color[i] = colormap[8]
        snowfall_color_z.append(dm_colormap[8])
        
    elif snow_z[i]>1.5 and snow_z[i]<=2.:
        
        snowfall_color_z.append(dm_colormap[9])
        
    elif snow_z[i]>2. and snow_z[i]<=2.5:
        
        snowfall_color_z.append(dm_colormap[10])
        
    elif snow_z[i]>2.5 and snow_z[i]<=3.:
        
        snowfall_color_z.append(dm_colormap[11])
        
    elif snow_z[i]>3. and snow_z[i]<=3.5:
        
        snowfall_color_z.append(dm_colormap[12])       
        
    elif snow_z[i]>3.5 and snow_z[i]<=4.:
        
        snowfall_color_z.append(dm_colormap[13])
        
    elif snow_z[i]>4. and snow_z[i]<=4.5:
        
        snowfall_color_z.append(dm_colormap[14]) 
        
        
        
#Now for KuPR        

snowfall_color_ku = []





for i in range(len(snow_ku)):

    if snow_ku[i]<=0.3:
        #z_color[i] = colormap[0]
        snowfall_color_ku.append(dm_colormap[0])
    
    elif snow_ku[i]>0.3 and snow_ku[i]<=0.4:
        #z_color[i] = colormap[1]
        snowfall_color_ku.append(dm_colormap[1])
    
    elif snow_ku[i]>0.4 and snow_ku[i]<=0.5:
        #z_color[i] = colormap[2]
        snowfall_color_ku.append(dm_colormap[2])
    
    elif snow_ku[i]>0.5 and snow_ku[i]<=0.6:
        #z_color[i] = colormap[3]
        snowfall_color_ku.append(dm_colormap[3])
    
    elif snow_ku[i]>0.6 and snow_ku[i]<=0.7:
        #z_color[i] = colormap[4]
        snowfall_color_ku.append(dm_colormap[4])
    
    elif snow_ku[i]>0.7 and snow_ku[i]<=0.8:
        #z_color[i] = colormap[5]
        snowfall_color_ku.append(dm_colormap[5])
        
    elif snow_ku[i]>0.8 and snow_ku[i]<=0.9:
        #z_color[i] = colormap[6]
        snowfall_color_ku.append(dm_colormap[6])
    
    elif snow_ku[i]>0.9 and snow_ku[i]<=1.:
        #z_color[i] = colormap[7]
        snowfall_color_ku.append(dm_colormap[7])
        
    elif snow_ku[i]>1. and snow_ku[i]<=1.5:
        #z_color[i] = colormap[8]
        snowfall_color_ku.append(dm_colormap[8])
        
    elif snow_ku[i]>1.5 and snow_ku[i]<=2.:
        
        snowfall_color_ku.append(dm_colormap[9])
        
    elif snow_ku[i]>2. and snow_ku[i]<=2.5:
        
        snowfall_color_ku.append(dm_colormap[10])
        
    elif snow_ku[i]>2.5 and snow_ku[i]<=3.:
        
        snowfall_color_ku.append(dm_colormap[11])
        
    elif snow_ku[i]>3. and snow_ku[i]<=3.5:
        
        snowfall_color_ku.append(dm_colormap[12])       
        
    elif snow_ku[i]>3.5 and snow_ku[i]<=4.:
        
        snowfall_color_ku.append(dm_colormap[13])
        
    elif snow_ku[i]>4. and snow_ku[i]<=4.5:
        
        snowfall_color_ku.append(dm_colormap[14]) 

        


 #Now for reflectivity bias:


bias_colormap = ['midnightblue', 'darkblue', 'mediumblue', 'blue', 'royalblue', 'cornflowerblue', 'lightcoral', 'coral', 'orangered','red','firebrick', 'darkred' ]

bias_color = []





for i in range(len(z_bias)):

    if z_bias[i]<=-5:
        #z_color[i] = colormap[0]
        bias_color.append(bias_colormap[0])
    
    elif z_bias[i]>-5 and z_bias[i]<=-4:
        #z_color[i] = colormap[1]
        bias_color.append(bias_colormap[1])
    
    elif z_bias[i]>-4 and z_bias[i]<=-3:
        #z_color[i] = colormap[2]
        bias_color.append(bias_colormap[2])
    
    elif z_bias[i]>-3 and z_bias[i]<=-2:
        #z_color[i] = colormap[3]
        bias_color.append(bias_colormap[3])
    
    elif z_bias[i]>-2 and z_bias[i]<=-1:
        #z_color[i] = colormap[4]
        bias_color.append(bias_colormap[4])
    
    elif z_bias[i]>-1 and z_bias[i]<=0:
        #z_color[i] = colormap[5]
        bias_color.append(bias_colormap[5])
        
    elif z_bias[i]>0 and z_bias[i]<=1:
        #z_color[i] = colormap[6]
        bias_color.append(bias_colormap[6])
    
    elif z_bias[i]>1 and z_bias[i]<=2:
        #z_color[i] = colormap[7]
        bias_color.append(bias_colormap[7])
        
    elif z_bias[i]>2. and z_bias[i]<=3:
        #z_color[i] = colormap[8]
        bias_color.append(bias_colormap[8])
        
    elif z_bias[i]>3. and z_bias[i]<=4.:
        
        bias_color.append(bias_colormap[9])
        
    elif z_bias[i]>4. and z_bias[i]<=5:
        
        bias_color.append(bias_colormap[10])
        
    elif z_bias[i]>=5 and z_bias[i]<=3.:
        
        bias_color.append(bias_colormap[11])
        
          
bias_colormap = ListedColormap(bias_colormap)        
        



#%%
import matplotlib

elev_angle = 0
label_size = 24
title_size = 28
tick_label_size = 20

#Setup plotting 
cmin = 0; cmax = 50; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_z,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = z_color, vmin = 0, vmax = 50, cmap = cmap, edgecolors = 'none')


for i in range(len(z_color)):
        if z_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX $0.5^{o}$ $Z_{H}$ 3/8/2018 ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label('dBZ', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)


#Now plot log10(nw)

cmin = 1; cmax = 5; cint = 0.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_z,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = nw_color, vmin = 1, vmax = 5, cmap = cmap, edgecolors = 'none')


for i in range(len(nw_color)):
        if nw_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX $0.5^{o}$ $log_{10}(N_{w})$ 3/8/2018 ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm $m^{-3}$]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)

###Now for Dm

cmin = 0.5; cmax = 3.5; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_dm,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(z_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = dm_color, vmin = 0.5, vmax = 3.5, cmap = cmap, edgecolors = 'none')


for i in range(len(dm_color)):
        if dm_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX $0.5^{o}$ $D_{m}$ 3/8/2018',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)




#Ground radar matchup snowfall rate


cmin = 0.; cmax = 4.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_dm,lut=nlevs)

colour_norm_object = matplotlib.colors.BoundaryNorm(boundaries = cmap_snow_array, ncolors  = 15)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(snowfall_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = snowfall_color, vmin = 0., vmax = 4., cmap = cmap_dm, edgecolors = 'none')


for i in range(len(snowfall_color)):
        if snowfall_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX S($Z_{H}$, $K_{DP}$) 3/8/2018',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm/hour]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)

#Now plot ice water content


cmin = 0.; cmax = 1.1; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_dm,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(iwc_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = iwc_color, vmin = 0., vmax = 1.1, cmap = cmap_dm, edgecolors = 'none')


for i in range(len(iwc_color)):
        if iwc_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX IWC($Z_{H}$, $K_{DP}$) 3/8/2018',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[$gm^{-3}$]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)


#%%

#Ground-based reflectivity

cmin = 0.; cmax = 4.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_dm,lut=nlevs)

colour_norm_object = matplotlib.colors.BoundaryNorm(boundaries = cmap_snow_array, ncolors  = 16)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(snowfall_color_z)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = snowfall_color_z, vmin = 0., vmax = 4., cmap = cmap_dm, edgecolors = 'none')


for i in range(len(snowfall_color_z)):
        if snowfall_color_z[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX S($Z_{H}$) 3/8/2018 ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm/hour]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)


# KuPR

cmin = 0.; cmax = 4.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=cmap_dm,lut=nlevs)

colour_norm_object = matplotlib.colors.BoundaryNorm(boundaries = cmap_snow_array, ncolors  = 16)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(snowfall_color_ku)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = snowfall_color_ku, vmin = 0., vmax = 4., cmap = cmap_dm, edgecolors = 'none')


for i in range(len(snowfall_color_ku)):
        if snowfall_color_ku[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX S(KuPR) 3/8/2018 ',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label(r'[mm/hour]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)


#%%


#Reflectivity Bias



cmin = -5.; cmax = 5.; cint = 1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name=bias_colormap,lut=nlevs)

colour_norm_object = matplotlib.colors.Normalize(vmin=cmin, vmax=cmax, clip=False)
scalar_mappable_object = plt.cm.ScalarMappable(cmap=cmap, norm=colour_norm_object)
scalar_mappable_object.set_array(bias_color)

figure_object, axes_object = plt.subplots(1, 1, figsize=(10, 10))

c = plt.scatter(longitude[:, elev_angle], latitude[:, elev_angle], c = bias_color, vmin = -5, vmax = 5, cmap = cmap, edgecolors = 'none')


for i in range(len(bias_color)):
        if bias_color[i] ==0: 
            c2 = plt.scatter(longitude[i,elev_angle], latitude[i,elev_angle], edgecolors = 'k', facecolors = 'none')
    
    
plt.xlabel('Longitude', size = label_size)
plt.ylabel('Latitude', size = label_size)
plt.xticks(size = tick_label_size)
plt.yticks(size = tick_label_size)
plt.title(r'KBOX $0.5^{o}$ Reflectivity Bias 3/2/2018',name='Calibri',size=26)

color_bar_object = plt.colorbar(ax=axes_object, mappable=scalar_mappable_object, orientation='vertical')

color_bar_object.set_label('[dBZ]', size = title_size)
color_bar_object.ax.tick_params(labelsize = label_size)

plt.show()
plt.close(figure_object)

