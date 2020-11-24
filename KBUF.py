#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 17:17:07 2020

@author: noahbrauer
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import netCDF4
from netCDF4 import Dataset
import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

import pyart
from mpl_toolkits.basemap import Basemap

file_zh = 'KBUF_N0R_20180302_040200_Zh.nc'
file_zdr = 'KBUF_N0X_20180302_040200_Zdr.nc'

nc_zh = Dataset(file_zh, 'r')
nc_zdr = Dataset(file_zdr,'r')

lat = nc_zh.variables['lat'][:]
lon = nc_zh.variables['lon'][:]

lat2,lon2 = np.meshgrid(lat,lon)
time = nc_zh.variables['time'][:]
z = nc_zh.variables['bref'][:]
zdr = nc_zdr.variables['zdr'][:]





def diff_reflect():
    diff_reflect_cdict ={
                'red':((0.000, 0.000, 0.000),
                           (0.333, 1.000, 1.000),
        				 (0.417, 0.000, 0.000),
                           (0.500, 0.000, 0.000),
                           (0.583, 1.000, 1.000),
                           (0.750, 1.000, 1.000),
                           (0.833, 1.000, 1.000),
                           (1.000, 1.000, 1.000)),
        		'green':	((0.000, 0.000, 0.000),
                           (0.333, 1.000, 1.000),
        				 (0.417, 0.000, 0.000),
                           (0.500, 1.000, 1.000),
                           (0.583, 1.000, 1.000),
                           (0.750, 0.000, 0.000),
                           (0.833, 0.000, 0.000),
                           (1.000, 1.000, 1.000)),
                  'blue': ((0.000, 0.000, 0.000),
                           (0.333, 1.000, 1.000),
        				 (0.417, 1.000, 1.000),
                           (0.500, 1.000, 1.000),
                           (0.583, 0.000, 0.000),
                           (0.750, 0.000, 0.000),
                           (0.833, 1.000, 1.000),
                           (1.000, 1.000, 1.000))}
    diff_reflect_coltbl = LinearSegmentedColormap('DIFF_REFLECT_COLTBL',diff_reflect_cdict)
    return diff_reflect_coltbl
zdrcolor = diff_reflect()




plt.figure(figsize=(14,14))

min_value = -5
max_value = 75
value_interval = 2.5
title_font_size = 28
label_size = 20

cmin = min_value; cmax = max_value; cint = value_interval; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-81,-76]); ylim = np.array([41.5,45])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(); m.drawcountries()

cs = m.contourf(lon2,lat2,z[0,:,:].T, clevs, cmap = 'pyart_NWSRef', extend = 'both')

m.drawcounties()

cbar = plt.colorbar(fraction=0.03)
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.title(r'KBUF 3/2 0402 UTC $Z_{H}$', size = title_font_size)
plt.show()



plt.figure(figsize=(14,14))

min_value = -2
max_value = 5
value_interval = 0.25
title_font_size = 28
label_size = 20

cmin = min_value; cmax = max_value; cint = value_interval; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-81,-76]); ylim = np.array([41.5,45])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(); m.drawcountries()

cs = m.contourf(lon2,lat2,zdr[0,:,:].T, clevs, cmap = zdrcolor, extend = 'both')

m.drawcounties()

cbar = plt.colorbar(fraction=0.03)
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dB]',size = label_size)
plt.title(r'KBUF 3/2 0402 UTC $Z_{DR}$', size = title_font_size)
plt.show()



