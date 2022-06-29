# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 17:50:52 2022

@author: yanan
"""
import pickle
import json
import pandas as pd

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)
    

filename='E:/NALMA_LIS/NALMA_LIS_matches.pkl'

M=load_obj(filename)

for i, m in M.items():
    m['LMA'].t=m['LMA'].t.astype('str') # JSON DOEST NOT ACCEPT pd.Timestamp,so we need to make it a string
    m['LMA']=m['LMA'].to_dict('list') # make it dict and each item in the dict is a list so that we could have one key for one list (array) correspondence 
    
    
    if 'LIS_events' in m.keys():
        #m['LIS_events']=m['LIS_events'][['time','lat','lon','radiance','id','_orig_id','parent_id']] # narrow down to fields that you actually needs
        m['LIS_events'].time=m['LIS_events'].time.astype('str')
        m['LIS_events']=m['LIS_events'].to_dict('list')
        
    if 'RS' in m.keys():
        m['RS'].t=m['RS'].t.astype('str')
        m['RS']=m['RS'].to_dict('list')
    
    if 'LIS_events_polygon_xy' in m.keys():    
        a=m['LIS_events_polygon_xy']
        b=[x.tolist() for x in a]
        m['LIS_events_polygon_xy']=b
        
    m['centroid pxpy']=m['centroid pxpy'].tolist()    
    m['lma_flash_centroid_location']=m['lma_flash_centroid_location'].tolist()


with open('NALMA_LIS.json', 'w') as fp:
    json.dump(
        M, 
        fp, 
        indent=4, 
        sort_keys=True
    )