# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 15:45:56 2022
@author: yanan

This script serves to find matched-passover LIS files that have optical lightning
that were detected from space over a lightning mapping array detection domain.

Given a folder with all LIS files that you want to perform a search, 
it saves the files names in CSV format. This script were kept as a separate one simply because the searching took a while. 
It would be more efficient if we just run such search once and save the matched LIS files for future usage.

In this example, we used NALMA and only search for LIS events within 120 km of NALMA center

"""

import os
from pyltg.core.lis import LIS
import pymap3d as pm
import numpy as np
import matplotlib.pyplot as plt



def haversine_distance(latlon1, latlon2):                                                                       
    """                                                                                                                 
    Calculates the distances between the points          
    They are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
    the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: N*2 nparray                                                                         
        latitudes:  N*2 nparray                                                                                                                                                
 
    Returns:                                                                                                            
        distance in km                                                                             
 
    """
    assert len(latlon1)==len(latlon2)    
                                                                                                             
    # if it is 1d array, make it 1*2
    if len(latlon1.shape) == 1:
        latlon1=latlon1.reshape(1,-1)
    if len(latlon2.shape) == 1:                                                                                       
        latlon2=latlon2.reshape(1,-1)                                                                        
 
    lon1 = latlon1[:,1]                                                                                      
    lat1 = latlon1[:,0]                                                                                       
    lon2 = latlon2[:,1]                                                                                        
    lat2 = latlon2[:,0]                                                                                         
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    km_array = 6371 * 2 * np.arcsin(np.sqrt(a))                                                                         
                                                                   
    return km_array



NALMA_coordinates=np.array([[34.8092586,  -87.0357225],
                            [34.6433808,  -86.7714025],
                            [34.7253536,  -86.6449781],
                            [34.6656331,  -86.3586129],
                            [34.7455622,  -86.5126506],
                            [34.9836496,  -86.8393545],
                            [34.8996936,  -86.5578487],
                            [34.6121906,  -86.5196873],
                            [34.5231382,  -86.9681644],
                            [34.6578669,  -87.3436469],
                            [35.1532036,  -87.0611744],
                            [35.0684567,  -86.5624089]])

#LMA_center=np.array([34.8,-86.85]) # NALMA
LMA_center=np.array([35.3,-98.5]) # OKLMA
distance_threshold=150


data_dir='E:/LIS_data/'
LMA_network_name='OKLMA'
fname_list=[]

for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if name[-3:]=='.nc':
            name=os.path.join(path, name).replace("\\" , "/")
            fname_list.append(name)
    
#keep the files with flashes within 80 km of LMA center
passover_fname_list=[]
num_flashes_within_lma=[]
num_null_files=0


matched_filenames=open("LIS_"+LMA_network_name+"_matched_filenames.txt", "w")

for i, fname in enumerate(fname_list):
    print(str(i+1)+'/'+str(len(fname_list)),fname)
    l = LIS(fname)
    
    # if this file has no events data, l.flashes will be empty list, and we will skip it
    if len(l.flashes)==0:
        num_null_files+=1
        continue
    
    f=l.flashes._data
    f_lat=f['lat'].values
    f_lon=f['lon'].values
    f_latlon=np.hstack((f_lat.reshape(-1,1),f_lon.reshape(-1,1)))
    
    n_flashes=len(f_lat)
    
    NALMA_center_duplicates=np.tile(LMA_center, (n_flashes, 1))
    f_d_2_lma=haversine_distance(f_latlon,NALMA_center_duplicates)
    
    num_f_within_lma=len(np.where(f_d_2_lma<distance_threshold)[0])
    
    if num_f_within_lma>0:
       passover_fname_list.append(fname)
       num_flashes_within_lma.append(num_f_within_lma)
       matched_filenames.write(fname + "\n")

matched_filenames.close()


    
    
    
    
    
    
    