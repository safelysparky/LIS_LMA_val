# -*- coding: utf-8 -*-
"""
Created on Tue May  3 11:02:34 2022

@author: yanan
"""



from pyltg.core.lis import LIS
import pyltg.core.lis as ltgLIS

import os
import pymap3d as pm
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import gzip
import datetime

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry import LinearRing
from scipy.spatial import ConvexHull
from matplotlib import cm
from astral.sun import sun
from astral import LocationInfo

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def find_time_ranges_of_LIS_events_over_LMA (l,LMA_center,distance_thres_km):
    E=l.events.data
    
    # Find the rought time range when LIS detected events within X km of the LMA center
    # Later this time range will be used to extract LMA and ENTLN data
    E_lat=E['lat'].values
    E_lon=E['lon'].values
    E_latlon=np.hstack((E_lat.reshape(-1,1),E_lon.reshape(-1,1)))
        
    E_d_2_lma=haversine_distance(E_latlon,LMA_center)
    
    # if LIS flash is within 80 km of LMA center
    close_E_idx=np.where(E_d_2_lma<distance_thres_km)[0]
    
    num_E_within_lma=len(close_E_idx)
    if num_E_within_lma>0:
    
        E_close=E.iloc[close_E_idx,:]
        
        #sort times stamps
        E_close_t=E_close['time'].sort_values()
        first_LIS_event_t=E_close_t.iloc[0]
        last_LIS_event_t=E_close_t.iloc[-1]
        print(first_LIS_event_t,last_LIS_event_t)
    else:
        first_LIS_event_t=[]
        last_LIS_event_t=[]
        print(f"no LIS events detected within {distance_thres_km} km of the LMA network ")
        
    return first_LIS_event_t,last_LIS_event_t

def find_sunrise_sunset_times(lat,lon,date_str):
    
    # create the sun observer location, the lat lon the centroid of the lma flash
    l = LocationInfo()
    l.latitude = lat
    l.longitude = lon
    
    now_date=datetime.datetime.strptime(date_str, '%Y%m%d').date() # the date of flash
    
    one_day=datetime.timedelta(days=1)
    yest_date=now_date-one_day
    tomor_date=now_date+one_day
    
    # get sunrise and sunset time for three consecutive days
    # simply to avoid any potential issue like sunrise and sunset might be on two different days
    # on top of this, in different time zone, sunset and sunrise order might change in UTC time 
    sunrise_t_list=[]
    sunset_t_list=[]
    for day in [yest_date,now_date,tomor_date]:
        s = sun(l.observer, date=day,tzinfo='UTC')
        # get the sunrise time and also remove tzinfo =UTC from the datetime object 
        # for easy comparison with f_t1 (with no tzinfo)
        sunrise_t=s['sunrise'].replace(tzinfo=None) 
        sunset_t=s['sunset'].replace(tzinfo=None)
        sunrise_t_list.append(sunrise_t)
        sunset_t_list.append(sunset_t)
        
        
    return sunrise_t_list,sunset_t_list
    
def determine_day_night_given_sun_times(sunrise_t_list,sunset_t_list,f_t1_tstamp_till_second):

    # f_t1_tstamp_till_us is the timing of first lma source in the flash, type is pd timestamp
    # need to convert it to datetime for comparison with sunrise and sunset
    f_t1_datetime=f_t1_tstamp_till_second.to_pydatetime()
    
    sun_rise_set_list=sunrise_t_list+sunset_t_list
    sun_rise_set_list.sort()
    
    for i, t in enumerate(sun_rise_set_list):
        if f_t1_datetime>=t and f_t1_datetime<=sun_rise_set_list[i+1]:
            break
        
    if (sun_rise_set_list[i] in sunrise_t_list)&(sun_rise_set_list[i+1] in sunset_t_list):
        dn="day"
        sunrise_diff_hours=np.around((f_t1_datetime-sun_rise_set_list[i]).total_seconds()/3600,1)
        sunset_diff_hours=np.around((f_t1_datetime-sun_rise_set_list[i+1]).total_seconds()/3600,1)
    
    if (sun_rise_set_list[i] in sunset_t_list)&(sun_rise_set_list[i+1] in sunrise_t_list):
        dn="night"
        sunset_diff_hours=np.around((f_t1_datetime-sun_rise_set_list[i]).total_seconds()/3600,1)
        sunrise_diff_hours=np.around((f_t1_datetime-sun_rise_set_list[i+1]).total_seconds()/3600,1)
        
    return sunrise_diff_hours, sunset_diff_hours, dn
    

def determine_flash_occur_in_day_or_night(f_t1_tstamp,ref_lat,ref_lon):
    
    date_str=''.join(str(f_t1_tstamp)[:10].split('-'))
    f_t1_tstamp_till_second=pd.Timestamp(str(f_t1_tstamp)[:19])
    
    sunrise_t_list,sunset_t_list=find_sunrise_sunset_times(ref_lat,ref_lon,date_str)
    
    sunrise_diff_hours, sunset_diff_hours, dn= determine_day_night_given_sun_times(sunrise_t_list,sunset_t_list,f_t1_tstamp_till_second)

    return sunrise_diff_hours, sunset_diff_hours, dn

            

from matplotlib.gridspec import GridSpec
def LIS_LMA_plot_layout():
    fig = plt.figure(figsize=(12,21))
    
    gs = GridSpec(42, 24)
    
    ax1 = plt.subplot(gs[0:8, :])
    ax2 = plt.subplot(gs[9:16, :])
    ax3 = plt.subplot(gs[17:, :])
    
    
    ax1.set_ylabel('h (m)')
    # ax1.set_xticks([])
    
    ax2.set_ylabel('Radiance')
    ax2.set_xlabel('Time (ms)')
    ax3.set_xlabel('EW (km)')
    ax3.set_ylabel('NS (km)')
    
    axs=[ax1,ax2,ax3]
    return fig,axs

def haversine_distance(latlon1, latlon2):                                                                       
    """                                                                                                                 
    Calculates the distances between the points          
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

def haversine_latlon_xy_conversion(latlon1, ref_latlon):                                                                       

    
    # if it is 1d array, make it 1*2
    if len(latlon1.shape) == 1:
        latlon1=latlon1.reshape(1,-1)
    if len(ref_latlon.shape) == 1:                                                                                       
        ref_latlon=ref_latlon.reshape(1,-1) 
    
    # lets duplicate the ref_latlon to make them equal length
    ref_latlon=np.tile(ref_latlon, (len(latlon1), 1))
    
    #convert lon to x 
    lon1 = latlon1[:,1]                                                                                      
    lat1 = ref_latlon[:,0]                                                                                       
    lon2 = ref_latlon[:,1]                                                                                        
    lat2 = ref_latlon[:,0]   
    
    #we need to cosider sign
    polarity=np.ones(lon1.shape)
    polarity[lon1<lon2]=-1
                                                                                      
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    x_km_array = 6371 * 2 * np.arcsin(np.sqrt(a)) 
    x_km_array = polarity*x_km_array                                                                        
    
    #convert lat to y 
    lon1 = ref_latlon[:,1]                                                                                     
    lat1 = latlon1[:,0]                                                                                      
    lon2 = ref_latlon[:,1]                                                                                        
    lat2 = ref_latlon[:,0]    

    #we need to cosider sign
    polarity=np.ones(lon1.shape)
    polarity[lat1<lat2]=-1
                                                                                     
 
    # Vectorized haversine calculation                                                                                  
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])                                                  
    a = np.sin((lat2 - lat1) / 2.0)**2 + (np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.0)**2)                 
    y_km_array = 6371 * 2 * np.arcsin(np.sqrt(a))  
    y_km_array = polarity*y_km_array 
    
    x_km_array=x_km_array.reshape(-1,1)
    y_km_array=y_km_array.reshape(-1,1)
    
    xy_km_array=np.hstack((x_km_array,y_km_array))
    
    return xy_km_array

def group_pd_tstamps_to_epoch(group_t_stamps):
    
    group_epoch=np.zeros(len(group_t_stamps),dtype=np.float64)
    for i, t_stamps in enumerate(group_t_stamps):
        t_epoch=t_stamps.value/(int(1e9))
        group_epoch[i]=t_epoch
        
    return group_epoch


def find_lma_data_line(txt_filename):
    # this function get the line number of "*** data ***" below which data are given
    # if it is gzipped file, we need to open it with gzip module
    if txt_filename[-2:] =='gz':
        with gzip.open(txt_filename) as f:
            for i, line in enumerate(f):            
                if b"*** data ***" in line: 
                    break
    else:
        # if it is not a gzipped file, regular .dat or txt file
        with open(txt_filename) as f:
            for i, line in enumerate(f):            
                if "*** data ***" in line: 
                    break  
    return i

def geodetic_to_enu(lla,ref_lat,ref_lon,ref_alt):
    
    # CAMMA 11 is the reference 
    # ref_lat=-31.6671004
    # ref_lon=-63.8828453
    # ref_alt=357.693
    ell_wgs84 = pm.Ellipsoid( "wgs84")
    
    if len(lla.shape)==1: # if xyz is only vector rather than an 2d array 
        lat=lla[0]
        lon=lla[1]
        alt=lla[2]
        
        x, y, z = pm.geodetic2enu(lat, lon, alt, ref_lat, ref_lon, ref_alt, ell=ell_wgs84, deg=True) 
        xyz=np.array([x,y,z])

    else:
        lat=lla[:,0]
        lon=lla[:,1]
        alt=lla[:,2]
    
        x, y, z = pm.geodetic2enu(lat, lon, alt, ref_lat, ref_lon, ref_alt, ell=ell_wgs84, deg=True) 
        x=x.reshape(-1,1)
        y=y.reshape(-1,1)
        z=z.reshape(-1,1)
        
        xyz=np.concatenate((x,y,z),axis=1)
        
    return xyz 


def read_lma_format_data_as_nparray_with_epoch_t(txt_filename,ref_lat,ref_lon,ref_alt):
    # eg., COLMA center
    # ref_lat=40.4463980
    # ref_lon=-104.6368130
    # ref_alt=1000
    fname_only=os.path.basename(txt_filename)
    lma_name, file_start_epoch_ns, duration_ns, date_str=LMA_fname_parse(fname_only)
    date_ref_t_ns=pd.Timestamp(date_str).value

    header_line=find_lma_data_line(txt_filename)+1
    # note here we only used the first six columns, because the sta msk in hex is not compatible
    L1=np.loadtxt(txt_filename,skiprows=header_line,usecols=(0,1,2,3,4,5))
    
    # For an file with no sources, loadtxt will return an empty umpy array
    # in this case, return S with empty list 
    if len(L1)==0:
        S=[]
        return S 
    
    # fix the one data line problem
    if len(L1.shape)==1:
        L1=L1.reshape(1,-1)
        
    lla=L1[:,1:4]

    rchi2=L1[:,4].reshape(-1,1)

    power=L1[:,5].reshape(-1,1)
    t=L1[:,0].reshape(-1,1)

    # convert tod to epoch t in seconds, note this only keep accuracy to us
    t=t+date_ref_t_ns/int(1e9)

    xyz=geodetic_to_enu(lla,ref_lat,ref_lon,ref_alt)
    xyz_km=xyz/1e3
    
    S=np.hstack((xyz_km,t,rchi2,power,lla))
    
    #!!!!!!!!!!delte source that are too high Or rchi2 too large!!!!!!!!!!!!
    altitude=lla[:,2]
    idx_keep=(altitude<=20*1e3)&(L1[:,4]<5)
    S=S[idx_keep,:]
    return S

def LMA_fname_parse(fname):
    #fname1: goesr_plt_NALMA_20170301_19_3600.dat.gz which only give hours without mmss
    #fname2: goesr_plt_COLMA_20170531_233000_0600.dat.gz which only give hhmmss
    fname_split=fname.split('_')
    for i, name_pieces in enumerate(fname_split):
        if 'LMA' in name_pieces:
            break
        if 'LYLOUT' in name_pieces: # some LMA files do not have LMA network name
            break
    lma_name=name_pieces
    date_str=fname_split[i+1]
    
    # note date str could be in form 20190801 or 190801 will have different values in pd.Timestamp()
    # to avoid this,add 20, here I boldly assume no data before 2000 will be used :-)
    if len(date_str)==6:
        date_str='20'+date_str
    
    hhmmss_str=fname_split[i+2]
    
    #sometimes, mm and ss are missing, in those case we need to pad zeros
    if len(hhmmss_str)==2:
        hhmmss_str=hhmmss_str+'0000'
    if len(hhmmss_str)==4:
        hhmmss_str=hhmmss_str+'00'
    
    
    file_start_epoch_ns=pd.Timestamp(date_str).value+(int(hhmmss_str[0:2])*3600+int(hhmmss_str[2:4])*60+int(hhmmss_str[4:]))*int(1e9)
    duration_str=fname_split[i+3].split('.')[0]    
    duration_ns=int(duration_str)*int(1e9) 
    
    return lma_name, file_start_epoch_ns, duration_ns, date_str

def find_involved_lma_file(NALMA_folder,LMA_NAME,first_LIS_event_t,last_LIS_event_t):
    
    first_t_ns=first_LIS_event_t.value
    last_t_ns=last_LIS_event_t.value
    
    involved_fnames=[]
    LMA_fnames=os.listdir(NALMA_folder)
    for fname in LMA_fnames:
        lma_name, file_start_epoch_ns, duration_ns, date_str=LMA_fname_parse(fname)
        #print(lma_name,date_str,file_start_tod,duration)
        file_end_epoch_ns=file_start_epoch_ns+duration_ns
        
        # first make sure LMA network name matches
        if (lma_name!=LMA_NAME) & (lma_name!="LYLOUT"):
            continue

        # flash comletely left of the lma file t range
        if (first_t_ns<file_start_epoch_ns) and (last_t_ns<file_start_epoch_ns):
            continue
        # flash comletely right of the lma file t range
        elif (first_t_ns>file_end_epoch_ns) and (last_t_ns>file_end_epoch_ns):
            continue
        else:
            involved_fnames.append(NALMA_folder+fname)
            
    return involved_fnames

def d_point_to_points_2d(event, events):
    dd=np.sqrt((events[:,0]-event[0])**2+(events[:,1]-event[1])**2)
    return dd

def extract_lma_data_in_LIS_passover_time_range(LMA_L1_folder,LMA_NAME,ref_lat,ref_lon,ref_alt,first_LIS_event_t,last_LIS_event_t):
    first_LIS_event_epoch=pd.Timestamp(first_LIS_event_t).value/(int(1e9))
    last_LIS_event_epoch=pd.Timestamp(last_LIS_event_t).value/(int(1e9))

    # find involved lma file names 
    involved_lma_files=find_involved_lma_file(LMA_L1_folder,LMA_NAME,first_LIS_event_t,last_LIS_event_t)
    
    if involved_lma_files==[]:
        print(f"No LMA files found for time range {first_LIS_event_t} - {last_LIS_event_t}")
        S_selected=[]
        return S_selected  # return here skip the rest of function
    
    #load the lma file and extract those only within the time ranges of passover LIS events:
    non_empty_file_num=0
    for j, lma_fname in enumerate(involved_lma_files):
        S_one_file=read_lma_format_data_as_nparray_with_epoch_t(lma_fname,ref_lat,ref_lon,ref_alt)
        
        if len(S_one_file)!=0:
            non_empty_file_num+=1
        else:
            print(f"The file {lma_fname} is empty!!!!!")
        
        if non_empty_file_num==1:
            S=S_one_file
        elif non_empty_file_num>1:
            # stack them if multiple files needs to be loaded
            S=np.vstack((S,S_one_file))
            
    # if all files are just empty        
    if non_empty_file_num==0:
        S_selected=[]
        return S_selected
            

    #shink it to the time range we needed, -/+1s of the LIS passover events time ranges:
    idx_keep=np.where((S[:,3]>first_LIS_event_epoch-1)&(S[:,3]<last_LIS_event_epoch+1))[0]
    
    if len(idx_keep)==0:
       S_selected=[]
       return S_selected
    
    S_selected=S[idx_keep,:]
    return S_selected


def flash_sorting(S,t_thres,d_thres,dd_method):
    #format of S [X,Y,Z,T,....] # make sure first 4 colums are xyzt.
    ## first sort S based on timings
    S=S[S[:,3].argsort()]
    flash_id_col=np.zeros(len(S))
    flash_id_col=flash_id_col.reshape(-1,1)
    S=np.hstack((S,flash_id_col))
    
    
    flash_info=np.zeros((len(S),4))
    # S=S[:50000]
    #sort sources into flashes
    for i, event in enumerate(S):
        if i==0:
            flash_id=0
            event[-1]=flash_id
            event=event.reshape(1,-1)
            # update buffer
            f_buffer=event
            # assign flash_id in the S
            S[0,-1]=flash_id
        else:
                    
            dt_last=event[3]-f_buffer[-1,3] # first event t minus the last t in the buffer
            if dt_last>t_thres: # if this source is more than t_thres away from all sources in the buffer, assign a new flash ID
                flash_id=flash_id+1
                event[-1]=flash_id
                f_buffer=np.vstack((f_buffer,event))
                
                #put flash_ID in the S
                S[i,-1]=flash_id
            
            else:
                dt=event[3]-f_buffer[:,3]
                if dd_method=='2d':
                    dd=d_point_to_points_2d(event,f_buffer)
                else:
                    dd=d_point_to_points_3d(event,f_buffer)
                
                bol=(dt<=t_thres)&(dd<=d_thres)
                
            
                if np.any(bol)==True: # if there is a source satisfy both conditions
                    sep_fac=(dt[bol]/0.01)+(dd[bol]/1) # separation factor, assume leader speed 10^5, 1 km equals 10 ms
                    match_s_idx=sep_fac.argmin() # idx of the matched cloest source
                    f_buffer_qua=f_buffer[bol]
                    match_s_fID=int(f_buffer_qua[match_s_idx,-1])
                    event[-1]=match_s_fID
                    
                    #put in the buffer
                    f_buffer=np.vstack((f_buffer,event))
                    
                    #put flash_ID in the S
                    S[i,-1]=match_s_fID
                    
                else:
                    # if no matched source found
                    flash_id=flash_id+1
                    #put flash_ID in the S
                    S[i,-1]=flash_id
                    # update buffer
                    event[-1]=flash_id
                    f_buffer=np.vstack((f_buffer,event))
                    
            ##############################################################################       
            # clean the buffer, clean buffer every 200 sources to make the function fast
            ##############################################################################                      
            if (i%200)==0:
            #if i>=0:
                f_buffer_t=f_buffer[:,3]
                df_t_buffer=event[3]-f_buffer_t
                old_events_idx=np.where(df_t_buffer>t_thres)[0]# those events that are too far away (temporally) away from the event
                kept_events_idx=np.where(df_t_buffer<=t_thres)[0]
                
                old_events_fID=np.unique(f_buffer[old_events_idx,-1])
                
                if len(old_events_fID)>0: # check if there are too old events
                    kept_events_fID=np.unique(f_buffer[kept_events_idx,-1])
                    
                    remove_tf=False
                    n_flash_2remove=0
                    for ID in old_events_fID:
                        if ID not in kept_events_fID: # if all events in the flash needs to be removed from the buffer, no events left in the buffer
                            n_flash_2remove=n_flash_2remove+1
                            remove_idx=np.where(f_buffer[:,-1]==ID)[0]
                            
                            if n_flash_2remove==1: # first flash found to be removed
                                remove_idx_total=remove_idx
                                remove_tf=True
                            else:
                                remove_idx_total=np.concatenate((remove_idx_total,remove_idx))
                                
                            ###############################################################    
                            ### calculate flash properties, flash start, end, n of sources
                            ###############################################################
                            flash_t1=f_buffer[remove_idx[0],3]
                            flash_t2=f_buffer[remove_idx[-1],3]
                            n_sources=len(remove_idx)
                            flash_info[int(ID),0:4]=np.array([flash_t1,flash_t2,n_sources,f_buffer[remove_idx[0],-1]])
                            
                    if remove_tf==True:
                        # print('delete')
                        f_buffer=np.delete(f_buffer, remove_idx_total, 0) 
                        
             # till the last event, calculte the flash properties of flashes left in buffer          
            if i==(len(S)-1):
                left_events_fID=np.unique(f_buffer[:,-1])
                for ID in left_events_fID:
                    
                    remove_idx=np.where(f_buffer[:,-1]==ID)[0]
                    
                    ###############################################################    
                    ### calculate flash properties, flash start, end, n of sources
                    ###############################################################
                    flash_t1=f_buffer[remove_idx[0],3]
                    flash_t2=f_buffer[remove_idx[-1],3]
                    n_sources=len(remove_idx)
                    flash_info[int(ID),0:4]=np.array([flash_t1,flash_t2,n_sources,f_buffer[remove_idx[0],-1]])
                 
    # delete those rows with zeros only, remove over allocated rows
    flash_info_sum=np.sum(flash_info,axis=1)
    flash_info=flash_info[flash_info_sum!=0]
    
    return S, flash_info

def filter_flashes_based_on_number_of_sources(S_sorted,flash_info,n_sources_thres):
    n_sources_col=flash_info[:,2]
    
    # flash info for big flash only
    big_flash_info=flash_info[n_sources_col>n_sources_thres,:]
    big_flash_fidx=big_flash_info[:,3]
    
    #here we keep S_sorted in big flashes only
    tf_row_keep=np.isin(S_sorted[:,-1],big_flash_fidx)
    S_sorted_big_flash=S_sorted[tf_row_keep,:]

    return S_sorted_big_flash, big_flash_info


def load_ENTLN_data(ENTLN_folder,ref_lat,ref_lon,ref_alt,first_LIS_event_t,last_LIS_event_t):
    e1_date_str=''.join(str(first_LIS_event_t)[:10].split('-'))
    e2_date_str=''.join(str(last_LIS_event_t)[:10].split('-'))

    ENTLN_file1=ENTLN_folder+e1_date_str+'.npy'
    if os.path.exists(ENTLN_file1) is True:
        EN=np.load(ENTLN_file1)
    else:
        print(f"EN data: {ENTLN_file1} does not exist, returned an empty array instead for cg_events")
        cg_events=np.array([])
        return cg_events
        

    # if LIS data span two days:
    if e2_date_str!=e1_date_str:
        ENTLN_file2=ENTLN_folder+e2_date_str+'.npy'
        if os.path.exists(ENTLN_file2) is True:
            EN_addtional=np.load(ENTLN_file2)
            EN=np.vstack((EN,EN_addtional))
        else:
            print(f"EN data: {ENTLN_file2} does not exist")

    #convert entln time to tod and lat lon to  xy in km
    EN_t_epoch=(EN[:,0]).reshape(-1,1)
    fake_EN_alt=np.zeros(len(EN_t_epoch)).reshape(-1,1) # assume all at zero 
    EN_lat=EN[:,1].reshape(-1,1)
    EN_lon=EN[:,2].reshape(-1,1)
    EN_type=EN[:,3].reshape(-1,1)
    EN_pc=(EN[:,4]/1e3).reshape(-1,1)
    
    EN_lla=np.hstack((EN_lat,EN_lon,fake_EN_alt)) 
    EN_xyz=geodetic_to_enu(EN_lla,ref_lat,ref_lon,ref_alt)
    EN_xy_km=EN_xyz[:,0:2]/1e3
    
    # EN_type are also for fake event ID, just to match format
    EN_events=np.hstack((EN_type,EN_xy_km,fake_EN_alt,EN_t_epoch,EN_pc,EN_lat,EN_lon)) 
    cg_events=EN_events[EN_type.flatten()==0,:]

    return cg_events     

def filter_flashes_within_FOV(S_sorted_big_flash, big_flash_info,l,M,in_FOV_lma_flash_no,EN_data_available,cg_events=[]):
    # by default, cg_events is empty, means no entln data were given
    
    for ii, info in enumerate(big_flash_info):
    
        flash_ID=info[3]
        sources_idx=np.where(S_sorted_big_flash[:,-1]==flash_ID)[0]

        # now the f is in the format of (x,y,z,t,chi,power,lla,flash_ID), in which xyz are in km
        f=S_sorted_big_flash[sources_idx]
                
        f_t=f[:,3]
        f_t1=np.min(f_t)
        f_t2=np.max(f_t)
        f_lat=f[:,6]
        f_lon=f[:,7]

        f_t1_tstamp=pd.to_datetime(f_t1, unit='s', origin='unix')
        f_t2_tstamp=pd.to_datetime(f_t2, unit='s', origin='unix')

        #get fov of LIS at flash t1 and flash t2
        fov1 = ltgLIS.get_fov(l.one_second, times=f_t1_tstamp)
        fov2 = ltgLIS.get_fov(l.one_second, times=f_t2_tstamp)
        
        #we need to make sure that LMA flashes are within the fovs
        fov1_lonlat=np.hstack((fov1[0].lon.reshape(-1,1),fov1[0].lat.reshape(-1,1)))
        fov2_lonlat=np.hstack((fov2[0].lon.reshape(-1,1),fov2[0].lat.reshape(-1,1)))
        
        #TODO: somehow some nan in the fov1 and fov2, skip this flash for now, need t double check that with PB
        if(np.isnan(fov1_lonlat[:,0]).any()) or (np.isnan(fov1_lonlat[:,1]).any()):
            print("The LIS fov Array contain NaN values")
            continue
        
        # make sure all LMA sources are within the fovs of LIS at flash t1 and t2
        in_fov1=check_if_within_polygon(fov1_lonlat,f_lat,f_lon)
        if in_fov1 is False:
            # print(ii,'out of FOV1')
            continue 
        
        in_fov2=check_if_within_polygon(fov2_lonlat,f_lat,f_lon)
        if in_fov2 is False:
            # print(ii,'out of FOV2')
            continue
                
        in_FOV_lma_flash_no+=1
        M[in_FOV_lma_flash_no]={}
        M[in_FOV_lma_flash_no]['LMA']=f
        M[in_FOV_lma_flash_no]['LMA_t1']=str(f_t1_tstamp)
        M[in_FOV_lma_flash_no]['LMA_t2']=str(f_t2_tstamp)
        
        # then look for cg events associated with the event, if ENTLN data is provided
        if EN_data_available==True: 

            flash_type=info[4] # flash number was given as a number
            if flash_type==1:
                flash_type_str='IC flash'
            elif flash_type==10:
                flash_type_str='Negative CG flash'
            elif flash_type==20:
                flash_type_str='Positive CG flash'
            else:
                flash_type_str='Bipolar flash'
            
            M[in_FOV_lma_flash_no]['type']=flash_type_str
            
            if flash_type!=1: # if the flash is a cg flash
                cg_idx=np.where(cg_events[:,-1]==ii)[0]
                cg_in_this_flash=cg_events[cg_idx]
                
                M[in_FOV_lma_flash_no]['RS']=cg_in_this_flash
                    
    return M,in_FOV_lma_flash_no

def get_lma_convex_hull_polygon(f_xy):
    
    # determine the convex hull of the LMA flash in 2D
    assert len(f_xy)>2 # MAKRE SURE lma SOURCES ARE MORE THAN 2
    
    lma_x=f_xy[:,0]
    lma_y=f_xy[:,1]
    hull = ConvexHull(f_xy)
    
    ## here we need to get the vertices of the convex hull and based on that create a polygon object
    hull_x=(lma_x[hull.vertices]).reshape(-1,1)
    hull_y=(lma_y[hull.vertices]).reshape(-1,1)
    hull_xy=np.hstack((hull_x,hull_y))
    rr=LinearRing(hull_xy)
    hull_polygon = Polygon(rr)

    return hull_polygon


def flash_type_assign(S_sorted_big_flash,cg_events,big_flash_info):
    
    # first check if cg_events is empty
    # if cg_events is empty, created a fake one so that all flashes will not match the cg and therefore all flashes are ICs
    if len(cg_events)==0:
        cg_events=np.array([-999999,-999999,-999999,-999999,-999999]).reshape(1,-1)
    
    
    # the column used to assign big_flash idx to a CG event, if the CG event is a match for a big flash
    # default for unmatches is -999
    parent_flash_idx_col=(np.ones(len(cg_events))*-999)

    # keep track of each cg_events separation factor relative to lma flash
    # because one CG event can be found as match for multiple lma flash
    # so we keep track of the sep facor and assign cg event to the lma flash with min sep_factor 
    cg_events_sep_factor=np.ones(len(cg_events))*99999999
    
    cg_t=cg_events[:,4]
    # big flash type default is 1, which is IC
    big_flash_type=np.ones(len(big_flash_info)).reshape(-1,1)
    
    # big flash no of strokes, default is 0
    big_flash_no_strokes=np.zeros(len(big_flash_info)).reshape(-1,1)
    
    # loop over each big flash to find the cloest cg in CG events data base.
    for ii,f_info in enumerate(big_flash_info):
        f_t0=f_info[0]-1*1e-3
        f_t1=f_info[1]+1*1e-3
        
        # first search based on time t
        idx_t_matched=np.where((cg_t-f_t0>=0)&(cg_t-f_t1<=0))[0]
        
        if len(idx_t_matched>0):

            #matched cg_events based on t
            matched_t_cg_event=cg_events[idx_t_matched,1:5]
            
            # get the sources for this flash
            tf_one_flash=(S_sorted_big_flash[:,-1]==big_flash_info[ii,3])
            S_one_flash=S_sorted_big_flash[tf_one_flash,:]
            
            # loop over each of the matched_t_cg_events, futher check spatial and temperal requirement
            for kk, one_cg_event in enumerate(matched_t_cg_event):
                row_no=idx_t_matched[kk] # the row num in cg_events
                
                dd=d_point_to_points_2d(one_cg_event,S_one_flash) # distance diff
                tt_diff=abs(S_one_flash[:,3]-one_cg_event[3]) # time diff in abs value
                
                
                bol=(dd<3)&(tt_diff<30*1e-3) # the parent flash has to have one sources within 30 ms and 3 km of the CG
                
                if np.any(bol)==True:
                    min_sep_factor=np.min(abs(tt_diff/0.01)+dd/1)  # separation factor, assume leader speed 10^5, 1 km equals 10 ms
                    if min_sep_factor<cg_events_sep_factor[row_no]:
                        cg_events_sep_factor[row_no]=min_sep_factor
                        parent_flash_idx_col[row_no]=ii # give the cg event the parent flash ID
    
    cg_events=np.hstack((cg_events,parent_flash_idx_col.reshape(-1,1))) 
    
    # now lets loop over big_flash_info again, after cg_events have been fianlly assigned to parent flashes
    for ii,f_info in enumerate(big_flash_info):
        
        matched_rows=np.where(parent_flash_idx_col==ii)[0]
        
        big_flash_no_strokes[ii]=len(matched_rows)
        
        if big_flash_no_strokes[ii]==0:
            big_flash_type[ii]=1 # np match found, IC
        else:
            pc=cg_events[matched_rows,5]
            n_pcg=np.sum(pc>0)
            n_ncg=np.sum(pc<0)
    
            if (n_pcg>0)&(n_ncg==0):
                big_flash_type[ii]=20
            elif (n_ncg>0)&(n_pcg==0):
                big_flash_type[ii]=10
            elif (n_pcg>0)&(n_ncg>0):
                big_flash_type[ii]=30
                
            
    ## attach the flash type and no of strokes to the big flash info
    big_flash_info=np.hstack((big_flash_info,big_flash_type,big_flash_no_strokes))
                    
    return big_flash_info,cg_events


def check_if_within_polygon(fov1_lonlat,f_lat,f_lon):
    polygon_fov1=Polygon(fov1_lonlat)
    in_polygon_tf=True
    for lat,lon in zip(f_lat,f_lon):
        point=Point(lon,lat)
        if polygon_fov1.contains(point) is False:
            in_polygon_tf=False
            break
            
    return in_polygon_tf

def convert_poly_latlon_to_xy(ev_poly,LMA_center):
    ev_poly_latlon=[]
    ev_poly_xy=[]
    for kk in range(len(ev_poly._paths)):
        poly_latlon=ev_poly._paths[kk].vertices[:,[1,0]]
        poly_xy=haversine_latlon_xy_conversion(poly_latlon,LMA_center)
        ev_poly_latlon.append(ev_poly._paths[kk].vertices[:,[1,0]])
        ev_poly_xy.append(poly_xy)
    
    return ev_poly_latlon,ev_poly_xy

def d_point_to_points_3d(event, events):
    dd=np.sqrt((events[:,0]-event[0])**2+(events[:,1]-event[1])**2+(events[:,2]-event[2])**2)
    return dd

def body_check(f,sep_thres):
    dd_col=np.zeros(len(f))
    for i,s in enumerate(f):
        dd=d_point_to_points_3d(s,f)
        dd=np.delete(dd,i) # delete the one with 0, distance to itself
        dd_col[i]=np.min(dd)
    
    tf_idx=dd_col<sep_thres
    f1=f[tf_idx,:] # delete those with large sepration to other sources
    
    return f1,tf_idx

def get_radiance_area_vs_time(E2):
    t_frame=np.array([])
    rad_frame=np.array([])
    area_frame=np.array([])
    
    for n in range(len(E2)):
        e_tstamp=E2.time.iloc[n]
        e_epoch=e_tstamp.value/(int(1e9))
        
        e_rad=E2.radiance.iloc[n]
        e_area=E2.footprint.iloc[n]
        
        frame_idx=np.where(t_frame==e_epoch)[0]
        
        if len(frame_idx)==0:
            t_frame=np.concatenate((t_frame,np.array([e_epoch])))
            rad_frame=np.concatenate((rad_frame,np.array([e_rad])))
            area_frame=np.concatenate((area_frame,np.array([e_area])))
        else:
            rad_frame[frame_idx]+=e_rad
            area_frame[frame_idx]+=e_area
            
    return t_frame,rad_frame,area_frame

def geolocate_all_pixels(t,one_second_df):
    #  Three element tuple given the interpolated position, velocity, and transformation matrix for the given time
    pol,vel,transform_matrix=ltgLIS.interp_one_sec(time=t,one_seconds=one_second_df)
    #  Here we geolocate all 128*128 pixels
    pxpy=np.array(np.meshgrid(np.arange(128), np.arange(128))).T.reshape(-1, 2)
    px_all=pxpy[:,0]
    py_all=pxpy[:,1]
    all_pixels_coordinates=ltgLIS.geolocate_pixels(pol,vel,transform_matrix,px_all,py_all)
    
    return all_pixels_coordinates,pxpy

def LIS_LMA_intesect(ev_poly_xy,hull_polygon_expanded):
    LIS_LMA_intesect_tf=False
    for poly in ev_poly_xy:
        poly_in_polygon_fmt=Polygon(poly)
        tf_pixel_intersect=poly_in_polygon_fmt.intersects(hull_polygon_expanded)
        
        if tf_pixel_intersect==True:
            LIS_LMA_intesect_tf=True
            break
        
    return LIS_LMA_intesect_tf

def expand_E2(E2,G,E):
    E2_parent_id=E2['parent_id'].values
    # prarent_group_indices = np.where(np.in1d(G_id, E2_parent_id))[0]
    E_parent_id=E['parent_id'].values
    E2_expanded_idx=np.where(np.in1d(E_parent_id, E2_parent_id))[0]
    E2_expanded=E.iloc[E2_expanded_idx]
    #replace E2 with E2 expanded
    E2=E2_expanded
    
    return E2

def LMA_LIS_match(M, batch_flashes_idx,l,LMA_center):
    for f_idx in batch_flashes_idx:
        # now the f is in the format of (x,y,z,t,chi,power,lla,flash_ID), in which xyz are in km
        f=M[f_idx]['LMA']
                
        f_t=f[:,3]
        f_t1=np.min(f_t)
        f_t2=np.max(f_t)
        f_lat=f[:,6]
        f_lon=f[:,7]
        f_latlon=np.hstack((f_lat.reshape(-1,1),f_lon.reshape(-1,1)))
        f_xy=haversine_latlon_xy_conversion(f_latlon,LMA_center)
        
        f_t1_tstamp=pd.to_datetime(f_t1, unit='s', origin='unix')
        # f_t2_tstamp=pd.to_datetime(f_t2, unit='s', origin='unix')
        
        # here we geolocate all pixels at the flash starting time, for later finding the hull centroid's pxpy
        all_pixels_coordinates,pxpy=geolocate_all_pixels(t=f_t1_tstamp,one_second_df=l.one_second)
        all_pixels_latlon=np.hstack((all_pixels_coordinates.lat.reshape(-1,1),all_pixels_coordinates.lon.reshape(-1,1)))
        all_pixels_xy=haversine_latlon_xy_conversion(all_pixels_latlon,LMA_center)
        
        # get the convexhull polygon of the lma flash
        # TODO: here need a line for hull in lat and lon
        hull_polygon=get_lma_convex_hull_polygon(f_xy)

        # expand the convex hull by 4 km to account of location offset of LIS
        hull_polygon_expanded = Polygon(hull_polygon.buffer(4.0).exterior) 
        
        #find the centroid of the LMA polygon (before expanding)
        hull_centroid_x,hull_centroid_y=hull_polygon.centroid.coords.xy
        
        #find the corresponding px and py for the centroid
        d_pixels_2_centroid=d_point_to_points_2d([hull_centroid_x,hull_centroid_y],all_pixels_xy)
        
        centroid_pixel_idx=d_pixels_2_centroid.argmin()
    
        centroid_pxpy=pxpy[centroid_pixel_idx]
        
        # detemine if it is day or night when this flash occurred       
        sunrise_diff_hours, sunset_diff_hours, dn = determine_flash_occur_in_day_or_night(f_t1_tstamp,LMA_center[0],LMA_center[1])
       
        M[f_idx]['flash area']=round(hull_polygon.area*100)/100
        M[f_idx]['centroid pxpy']=centroid_pxpy
        M[f_idx]['dn']= dn
        M[f_idx]['lma_flash_centroid_location']=np.array([hull_centroid_x,hull_centroid_y])
        
        # Seach for LIS events:
        # first lets narrow-down the events data based on temporal and spatial span of lma data        
        #temporal part
        G=l.groups.data
        E=l.events.data
        E_epoch=  group_pd_tstamps_to_epoch(E.time)
        idx_keep_temporal=np.where((E_epoch>=f_t1-2*1e-3)&(E_epoch<=f_t2+2*1e-3))[0]
        if len(idx_keep_temporal)==0:
            #print(f"base on time span of LMA flash no {f_idx}, no LIS events found")
            M[f_idx]['LIS detection']= False
            continue
        
        E1=E.iloc[idx_keep_temporal]
        
        #spatial part
        idx_keep_spatial=np.array([])
        for k in range(len(E1)):
            e_lat=E1.lat.iloc[k]
            e_lon=E1.lon.iloc[k]
            e_latlon=np.array([e_lat,e_lon])
            
            # the distance between an LIS event to all sources in a lma flash
            d_e_2_lma=haversine_distance(e_latlon, f_latlon)
            
            #we only keep events that are within 10 km of any lma sources
            if np.min(d_e_2_lma)<10:
                idx_keep_spatial=np.concatenate((idx_keep_spatial,np.array([k])))
                
        if len(idx_keep_spatial)==0:
            #print(f"base on spatial span of LMA flash no {f_idx}, no LIS events found")
            M[f_idx]['LIS detection']= False
            continue
        
        
        E2=E1.iloc[idx_keep_spatial]
        # optional, we could add addtional events that were not in E2 but in parent groups 
        E2=expand_E2(E2,G,E)
        
        # SORT E2 by time
        E2.sort_values(by=['time'])
        
        # get polygonS of all events in E2
        ev_poly_latlon,ev_poly_xy=get_polygons_of_events_in_E2(E2, l, LMA_center)
        
        # Now we can check if the expanded hull_polygon intersect with any pixels in E2
        LIS_LMA_intesect_tf= LIS_LMA_intesect(ev_poly_xy,hull_polygon_expanded)
        if LIS_LMA_intesect_tf ==True:
            # print('DETECTED')
            M[f_idx]['LIS detection']= True
            M[f_idx]['LIS_events']=E2
            M[f_idx]['LIS_events_polygon_xy']=ev_poly_xy
        else:
            # print('MISSED')
            M[f_idx]['LIS detection']= False
                    
    return M

def get_polygons_of_events_in_E2(E2,l,LMA_center):
    ev_poly_xy_E2=[]
    ev_poly_latlon_E2=[]
    
    last_e_stamp=pd.Timestamp('19700101')#just give a fake timestamp for comparison
    
    for n in range(len(E2)):
        e_tstamp=E2.time.iloc[n]
        e=pd.DataFrame([E2.iloc[n]])
        if e_tstamp!=last_e_stamp:# only run interp_one_sec when e_stamp changes, because it is slow
            pos, vel, xform_matrix = ltgLIS.interp_one_sec(e_tstamp, one_seconds=l.one_second)
            last_e_stamp=e_tstamp
        e_px = (e.px.values[:, np.newaxis] + np.array([-0.5, -0.5, 0.5, 0.5])).ravel()
        e_py = (e.py.values[:, np.newaxis] + np.array([-0.5, 0.5, 0.5, -0.5])).ravel()
        lla1 = ltgLIS.geolocate_pixels(pos, vel, xform_matrix, e_px, e_py)
        
        for key in lla1.__dict__.keys():
            setattr(lla1, key, getattr(lla1, key).reshape((-1, 4)))
        
        ev_poly1 = ltgLIS.event_poly(e, corners=lla1, latlon=True, fill=True)
        ev_poly_latlon1,ev_poly_xy1 = convert_poly_latlon_to_xy(ev_poly1,LMA_center)
        
        ev_poly_xy_E2.append(ev_poly_xy1[0])
        ev_poly_latlon_E2.append(ev_poly_latlon1[0])
        
    return ev_poly_latlon_E2,ev_poly_xy_E2



def make_arrays_into_dateframe(M):
    for key, m in M.items():
        lma_ary=m['LMA']
        lma_ary1=lma_ary[:,[0,1,8,3,4,5,6,7]]
        lma_ary1[:,2]=lma_ary1[:,2]/1e3 # convert lma altitude from m to km
        lma_df=pd.DataFrame(lma_ary1,columns=['x','y','altitude','t','chi','power','lat','lon'])
        lma_df.t=pd.to_datetime(lma_df.t,unit='s', origin='unix')
        m['LMA']=lma_df
        
        if 'RS' in m.keys():
            rs_ary=m['RS']
            rs_ary1=rs_ary[:,[1,2,3,4,5,6,7]]
            rs_df=pd.DataFrame(rs_ary1, columns = ['x','y','z','t','peak current','lat','lon'])
            rs_df.t=pd.to_datetime(rs_df.t, unit='s', origin='unix')
            m['RS']=rs_df
        
    return M


# def  main():

LMA_center=np.array([34.8,-86.85])
ref_lat=LMA_center[0]
ref_lon=LMA_center[1]
ref_alt=0


csv_fname ="C:/Users/yanan/Desktop/git_local/LIS_LMA_val/LIS_NALMA_matched_filenames.txt"

LMA_NAME='NALMA'
LMA_L1_folder='E:/NALMA_LIS/LMA data/'
ENTLN_folder='E:/NALMA_LIS/ENTLN data/'
fig_folder='C:/Users/yanan/Desktop/LIS_fig/'

EN_data_available=True

M={} # initialize the dictionary for LMA-LIS match
in_FOV_lma_flash_no=-1 # the initial value that will be used for LMA-LIS match dict key

df=pd.read_csv(csv_fname,names=['LIS_fnames','num_of_flashes','t1','t2'])

for i, row in df.iterrows():
    
    # if i+1!=1:
    #     continue
    fname=row['LIS_fnames']
    
    print(str(i+1)+'/'+str(len(df)),fname)

    ########################################################################
    # load the ~90 min LIS data
    ########################################################################
    l = LIS(fname)
    
    # find LIS events within X-km of LMA network center, the time ranges of those LIS events
    # will be used to extract LMA and ENTLN data later
    distance_thres_km=150
    first_LIS_event_t,last_LIS_event_t=find_time_ranges_of_LIS_events_over_LMA (l,LMA_center,distance_thres_km)
    # if no LIS events found, will return a empty list, then this file will be skipped
    if first_LIS_event_t==[]:
        continue 
    
    ##################################################################################
    # lets extract lma and entln data based on the passover LIS data time range:
    ##################################################################################
    # exrtact lma data in the time range:
    S_selected=extract_lma_data_in_LIS_passover_time_range(LMA_L1_folder,LMA_NAME,ref_lat,ref_lon,ref_alt,first_LIS_event_t,last_LIS_event_t)
    if len(S_selected)==0:
        continue
    
    #sort lma sources into flashes
    # flash sorting parameters:
    t_thres=0.2 # unit second
    d_thres=5 # unit km
    # S_sorted format: x(km),y(km),z(km),t(epoch in sec),rchi2,power,lla,flash_id
    # flash_info format: flash_t1,flash_t2,no_lma_sources,flash_id
    S_sorted,flash_info=flash_sorting(S_selected,t_thres,d_thres,dd_method='2d')

    # Here we keep only flashes with number of sources more than a threshold
    # since flashes with just a few sources could be just noises
    n_sources_thres=20
    S_sorted_big_flash, big_flash_info=filter_flashes_based_on_number_of_sources(S_sorted,flash_info,n_sources_thres)
    
    if len(S_sorted_big_flash)==0:
        print(f"No flashes more than {n_sources_thres} lma sources")
        continue
    
    # extract ENLTN data, Note include EN data is optional
    # note in this analysis, EN data is in numpy array format, and has already been filtered 
    # with temporal and spatial criteria via the data request, each EN data is saved as a file by date 
    if EN_data_available == True:
        cg_events=load_ENTLN_data(ENTLN_folder,ref_lat,ref_lon,ref_alt,first_LIS_event_t,last_LIS_event_t)

        # assign a flash as IC or CG based on if a ENTLN CG is present in the LMA flash:
        # a CG is considered belonging to the lma flash if the parent LMA flash have 
        # more than one sources within 30 ms and 3 km of the CG
        # big flash info format: flash_t0, flash_t1, num_sources, flash_idx, flash_type, number of strokes
        big_flash_info,cg_events=flash_type_assign(S_sorted_big_flash,cg_events,big_flash_info)
    
    
    #  further filter out LMA flashes that were out of the LIS fov and put those in FOV in a dictionary
    # in addtion, if entln cg data are provided, RS will be assgined to LMA flashes and flash type will be given
    M, in_FOV_lma_flash_no_updated =filter_flashes_within_FOV(S_sorted_big_flash, big_flash_info, l, M, in_FOV_lma_flash_no,EN_data_available, cg_events)
    
    # the flashes idx (key in M dictionary) of this batch LIS file
    batch_flashes_idx=np.arange(in_FOV_lma_flash_no+1,in_FOV_lma_flash_no_updated+1)
    
    # Seach LIS events for each of the LMA flash in FOV
    # If any LIS pixel intersect with the LMA flash convex hull 
    # then this flash is regarded as detected by LIS, and all LIS events will be saved
    # Some addtional properties of LMA flashes are determined
    # flash area, centroid of the convext hull as well as its pxpy
    # if a flash is in day or night 
    M=LMA_LIS_match(M, batch_flashes_idx,l,LMA_center)
    
    # update this number so that M can append
    in_FOV_lma_flash_no=in_FOV_lma_flash_no_updated
    


# Option to make lma array and entln array to be pandas df to increase the readibility
# and also remove some temp columns that no longer needed
M=make_arrays_into_dateframe(M)



pkl_name='E:/NALMA_LIS/'+LMA_NAME+'_LIS_matches.pkl'
save_obj(M,pkl_name)





 