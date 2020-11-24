#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 17:46:19 2020

@author: noahbrauer
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import Dataset
import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

import pyart





file = 'KBUF20180302_045130_V06'

radar = pyart.io.read(file)


    
xsect = pyart.util.cross_section_ppi(radar,[180])






display = pyart.graph.RadarDisplay(xsect)
fig = plt.figure(figsize = (10,10))
display.plot('differential_reflectivity', 0, vmin=-1, vmax=3., colorbar_label = 'dB')
plt.xlim(0,100)
plt.ylim(0,10)
plt.tight_layout()
plt.title(r'KBUF 3/2/2017 0402 UTC $Z_{DR}$', size = 20)
plt.ylabel('Altitude (km)',fontsize = 20)
plt.xlabel('Distance from radar (km)',fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.show()

