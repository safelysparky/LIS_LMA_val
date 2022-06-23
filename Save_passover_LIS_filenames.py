# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 15:45:56 2022
@author: yanan

This script serves to find matched-passover LIS files that have optical lightning
that were detected from space over a lightning mapping array detection domain.

Given a folder with all LIS files that you want to perform a search, 
it saves the files names in CSV format. 

This script were kept as a separate one from the LIS_LMA_match.py simply because the searching took a while for years of LIS files 
It would be more efficient if we just run such search once and save the matched LIS filenames for future usage.

In this example, we used NALMA and only search for LIS events within 150 km of NALMA center

"""

import os
from pyltg.core.lis import LIS
import numpy as np
import pandas as pd


def haversine_distance(latlon1, latlon2):                                                                       
    """                                                                                                                 
    Calculates the distances between points          
    They are calculated using a vectorized haversine calculation the great circle distance between two arrays of points on   
    the earth (specified in decimal degrees). All args must be of equal length.                                         
 
    Args:                                                                                                               
        longitudes: N*2 or 1*2 nparray                                                                         
        latitudes:  N*2 or 1*2 nparray                                                                                                                                                
 
    Returns:                                                                                                            
        distance in km                                                                             
 
    """
    
    # if it is 1d array, make it 1*2
    if len(latlon1.shape) == 1:
        latlon1=latlon1.reshape(1,-1)
    if len(latlon2.shape) == 1:                                                                                       
        latlon2=latlon2.reshape(1,-1) 
    
    # then if one is 1*2 but the other is N*2, lets duplicate the short one to make them equal length
    if (len(latlon1)==1)&(len(latlon2)>1):
        latlon1=np.tile(latlon1, (len(latlon2), 1))
    elif (len(latlon2)==1)&(len(latlon1)>1):
        latlon2=np.tile(latlon2, (len(latlon1), 1))
        
        
    
    assert len(latlon1)==len(latlon2)    
                                                                                                             
                                                                       
 
    lon1 = latlon1[:,1]                                                                                      
    lat1 = latlon1[:,0]                                                                                       
    lon2 = latlon2[:,1]                                                                                        
    lat2 = latlon2[:,0]                                                                                         
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    km_array = 6371 * 2 * np.arcsin(np.sqrt(a))                                                                         
                                                                   
    return km_array



LMA_center=np.array([34.8,-86.85]) # NALMA
distance_threshold=150 # Distance threshold set for looking for LIS events


data_dir='E:/LIS_data/'
LMA_name='NALMA'


fname_list=[]

#grab .nc files' full path in the folder, assume each filename has "LIS" in it:
for path, subdirs, files in os.walk(data_dir):
    for name in files:
        if (name[-3:]=='.nc') & ('LIS' in name) :
            name=os.path.join(path, name).replace("\\" , "/")
            fname_list.append(name)
    
#keep the files with flashes within 80 km of LMA center
passover_fname_list=[]
num_flashes_within_lma=[]
num_null_files=0


matched_filenames=open("LIS_"+LMA_name+"_matched_filenames.txt", "w")



for i, fname in enumerate(fname_list):
    print(str(i+1)+'/'+str(len(fname_list)),fname)
    l = LIS(fname)
    
    # if this file has no events data, l.flashes will be empty list, and we will skip it
    if len(l.flashes.data)==0:
        num_null_files+=1
        continue
    
    f=l.flashes.data
    f_lat=f['lat'].values
    f_lon=f['lon'].values
    f_latlon=np.hstack((f_lat.reshape(-1,1),f_lon.reshape(-1,1)))
    
    n_flashes=len(f_lat)
    
    f_d_2_lma=haversine_distance(f_latlon,LMA_center)
    within_thres_idx=np.where(f_d_2_lma<distance_threshold)[0]
    
    
    num_f_within_lma=len(within_thres_idx)
    
    if num_f_within_lma>0:
        
        t_close_flashes=f['time'].iloc[within_thres_idx]
        t1_close_flashes=t_close_flashes.min()
        t2_close_flashes=t_close_flashes.max()
        
        passover_fname_list.append(fname)
        num_flashes_within_lma.append(num_f_within_lma)
        
        # each column: fname full path, num_LIS_flashes, first_LIS_flash_time, last_LIS_flash_time
        matched_filenames.write(f"{fname},{num_f_within_lma},{t1_close_flashes},{t2_close_flashes}\n")
        

matched_filenames.close()

    
    
    
    
    
    
    