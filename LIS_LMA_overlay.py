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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import gzip
from pyproj import Transformer
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry import LinearRing
from scipy.spatial import ConvexHull
from matplotlib import cm


lonlat_to_webmercator = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
def latlon_to_Mercator(lon, lat):
    x, y = lonlat_to_webmercator.transform(lon, lat)
    return x, y


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


from matplotlib.gridspec import GridSpec
def LIS_plot_layout():
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

def group_t_stamps_to_tod(group_t_stamps):
    first_group_t_stamp=str(group_t_stamps.iloc[0])[:10] #eg '2017-03-01'
    ref_unix_ns=pd.Timestamp(first_group_t_stamp).value # use it as reference
    group_tod=[]
    for t_stamps in group_t_stamps:
        t_tod=(t_stamps.value-ref_unix_ns)/1e9
        group_tod.append(t_tod)
        
    return group_tod

def LIS_tstamp_to_tod(tstamp):
    tstamp_str=str(tstamp)[:10] #eg '2017-03-01'
    ref_unix_ns=pd.Timestamp(tstamp_str).value # use it as reference
    t_tod=(tstamp.value-ref_unix_ns)/1e9
    
    return t_tod

def tod_2_tstamp(tod):
    h=int(tod/3600)
    m=int((tod-h*3600)/60)
    s=tod-h*3600-m*60
    
    if h<10:
        h_str='0'+str(h)
    else:
        h_str=str(h)
    
    if m<10:
        m_str='0'+str(m)
    else:
        m_str=str(m)
        
    if s<10:
        s_str='0'+str(s)
    else:
        s_str=str(s)

    t_stamp=h_str+':'+m_str+':'+s_str
    
    return t_stamp

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


def read_lma_format_data_as_nparray(txt_filename,ref_lat,ref_lon,ref_alt):
    #COLMA center
    # ref_lat=40.4463980
    # ref_lon=-104.6368130
    # ref_alt=1000
    header_line=find_lma_data_line(txt_filename)+1
    # note here we only used the first six columns, because the sta msk in hex is not compatible
    L1=np.loadtxt(txt_filename,skiprows=header_line,usecols=(0,1,2,3,4,5))
    
    # fix the one data line problem
    if len(L1.shape)==1:
        L1=L1.reshape(1,-1)
        
    lla=L1[:,1:4]
    
    rchi2=L1[:,4].reshape(-1,1)
    power=L1[:,5].reshape(-1,1)
    t=L1[:,0].reshape(-1,1)
    xyz=geodetic_to_enu(lla,ref_lat,ref_lon,ref_alt)
    xyz_km=xyz/1e3
    
    S=np.hstack((xyz_km,t,rchi2,power,lla))
    
    #!!!!!!!!!!delte source that are too high Or rchi2 too large!!!!!!!!!!!!
    altitude=lla[:,2]
    idx_keep=(altitude<=20*1e3)&(L1[:,4]<5)
    S=S[idx_keep,:]
    return S

def LMA_fnames_parse(fname):
    #fname1: goesr_plt_NALMA_20170301_19_3600.dat.gz which only give hours without mmss
    #fname2: goesr_plt_COLMA_20170531_233000_0600.dat.gz which only give hhmmss
    fname_split=fname.split('_')
    for i, name_pieces in enumerate(fname_split):
        if 'LMA' in name_pieces:
            break
    lma_name=name_pieces
    date_str=fname_split[i+1]
    hhmmss_str=fname_split[i+2]
    
    #sometimes, mm and ss are missing, in those case we need to pad zeros
    if len(hhmmss_str)<=6:
        if len(hhmmss_str)==2:
            hhmmss_str=hhmmss_str+'0000'
        if len(hhmmss_str)==4:
            hhmmss_str=hhmmss_str+'00'
    
    
    file_start_tod=int(hhmmss_str[0:2])*3600+int(hhmmss_str[2:4])*60+int(hhmmss_str[4:])
    duration_str=fname_split[i+3].split('.')[0]    
    duration=int(duration_str)  
    
    return lma_name, date_str,file_start_tod, duration

def find_involved_lma_file(NALMA_folder,LMA_NAME,DATE_STR,flash_t0,flash_t1):

    involved_fnames=[]
    LMA_fnames=os.listdir(NALMA_folder)
    for fname in LMA_fnames:
        lma_name, date_str,file_start_tod, duration=LMA_fnames_parse(fname)
        print(lma_name,date_str,file_start_tod,duration)
        file_end_tod=file_start_tod+duration
        
        # first make sure LMA network name matches
        if LMA_NAME != lma_name:
            continue
        
        # then we mak esure date matches
        if DATE_STR != date_str:
            continue
        
        # case1 flash contained in one LMA file
        if (flash_t0>file_start_tod)&(flash_t0<file_end_tod):
            involved_fnames.append(NALMA_folder+fname)
            break
        
        # case2 flash intersecrted with LMA file beginning
        if (flash_t0<file_start_tod)&(flash_t1>file_start_tod):
            involved_fnames.append(NALMA_folder+fname)
         
        # case3 flash intersecrted with LMA file ending
        if (flash_t0<file_end_tod)&(flash_t1>file_end_tod):
            involved_fnames.append(NALMA_folder+fname)
        
    return involved_fnames

def d_point_to_points_2d(event, events):
    dd=np.sqrt((events[:,0]-event[0])**2+(events[:,1]-event[1])**2)
    return dd


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
                    
            ####################################################       
            # clean the buffer, clean buffer every 200 sources
            ####################################################                      
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
                 
    # delete those row with zeros only
    flash_info_sum=np.sum(flash_info,axis=1)
    flash_info=flash_info[flash_info_sum!=0]
    
    return S, flash_info

def flash_type_assign(S_sorted_big_flash,cg_events,big_flash_info):
    
    # the column used to assign big_flash idx to a CG event, if the CG event is a match for a big flash
    # default for unmatche is -999
    parent_flash_idx_col=(np.ones(len(cg_events))*-999)
    # cg_events=np.hstack((cg_events,big_flash_idx_col))
    cg_events_sep_factor=np.ones(len(cg_events))*99999999
    
    cg_t=cg_events[:,4]
    # big flash type default is 1, which is IC
    big_flash_type=np.ones(len(big_flash_info)).reshape(-1,1)
    
    # big flash no of strokes, default is 0
    big_flash_no_strokes=np.zeros(len(big_flash_info)).reshape(-1,1)
    
    # loop over each big flash
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
            
            # import pdb;pdb.set_trace()
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
    
    # now lets loop over big_flash_info again, after cg_events have been assigned to parent flashes
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

def convert_poly_latlon_to_xy(ev_poly):
    ev_poly_latlon=[]
    ev_poly_xy=[]
    for kk in range(len(ev_poly._paths)):
        poly_latlon=ev_poly._paths[kk].vertices[:,[1,0]]
        poly_xy=haversine_latlon_xy_conversion(poly_latlon,NALMA_center)
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
        e_tod=LIS_tstamp_to_tod(e_tstamp)
        
        e_rad=E2.radiance.iloc[n]
        e_area=E2.footprint.iloc[n]
        
        frame_idx=np.where(t_frame==e_tod)[0]
        
        if len(frame_idx)==0:
            t_frame=np.concatenate((t_frame,np.array([e_tod])))
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

NALMA_center=np.array([34.8,-86.85])
ref_lat=NALMA_center[0]
ref_lon=NALMA_center[1]
ref_alt=0


txt_file = open("C:/Users/yanan/Desktop/git_local/LIS_LMA_val/LIS_NALMA_matched_filenames.txt", "r")

fname_list=txt_file.read().splitlines() 

LMA_NAME='NALMA'
NALMA_folder='E:/NALMA_LIS_project/'
ENTLN_folder='E:/ENTLN_LIS_project/'


fig_folder='C:/Users/yanan/Desktop/LIS_fig/'


loaded_lma_fname=''

# We're going to be doing some plotting, so set up some variables
map_proj = ccrs.Mercator()
ll_proj_geo = ccrs.Geodetic()
ll_proj = ccrs.PlateCarree()

# flash sorting parameters:
t_thres=0.2 # unit second
d_thres=5 # unit km
dd_method='2d'

#get the colormap
cmap=cm.get_cmap('binary',100)

# each row correspond to one LMA flash in FOV OF lis, lma_flash_no, Detected by LIS (true or false), X, Y
lma_flash_no=0
DE_x_y=[]


for i, fname in enumerate(fname_list[0:1]):
    # print(str(i+1)+'/'+str(len(fname_list)),fname)
    l = LIS(fname)
    
    
    F=l.flashes._data
    G=l.groups.data
    E=l.events.data
    
    F_lat=F['lat'].values
    F_lon=F['lon'].values
    F_latlon=np.hstack((F_lat.reshape(-1,1),F_lon.reshape(-1,1)))
    
    n_flashes=len(F_lat)
    
    F_d_2_lma=haversine_distance(F_latlon,NALMA_center)
    
    # if LIS flash is within 80 km of LMA center
    close_F_idx=np.where(F_d_2_lma<80)[0]
    num_F_within_lma=len(close_F_idx)
    
    F_close=F.iloc[close_F_idx,:]
    
    #sort times stamps
    F_close_t=F_close['time'].sort_values()
    first_close_flash_t=F_close_t.iloc[0]
    last_close_flash_t=F_close_t.iloc[-1]
    print(first_close_flash_t,last_close_flash_t)
    print('_____________________________________________________________')
    
    ##################################################################################
    #now we have f_close, lets extract lma and entln data in the f_close time ranges:
    ##################################################################################
    
    # lets first extact lma data:
    
    # find the involved lma files
    flash_DATE_STR=''.join(str(first_close_flash_t)[:10].split('-'))
    
    ref_t_ns_4_tod=pd.Timestamp(flash_DATE_STR).value # use it as reference when calculate tod
    
    first_close_flash_tod=LIS_tstamp_to_tod(first_close_flash_t)
    last_close_flash_tod=LIS_tstamp_to_tod(last_close_flash_t)
    involved_lma_files=find_involved_lma_file(NALMA_folder,LMA_NAME,flash_DATE_STR,first_close_flash_tod,last_close_flash_tod)
    
    #load the lma file and extract those only within the f_close time ranges:
    for j, lma_fname in enumerate(involved_lma_files):
        S_one_file=read_lma_format_data_as_nparray(lma_fname,ref_lat,ref_lon,ref_alt)
        
        if j==0:
            S=S_one_file
        else:
            # stack them if multiple files needs to be loaded
            S=np.vstack((S,S_one_file))
    
    #shink it to the time range we needed, -/+3s of the f_close ranges:
    idx_keep=np.where((S[:,3]>first_close_flash_tod-3)&(S[:,3]<last_close_flash_tod+3))[0]
    S_selected=S[idx_keep,:]
    
    #now lets sort lma sources into flashes
    S_sorted,flash_info=flash_sorting(S_selected,t_thres,d_thres,dd_method='2d')
    

    # lets then extract ENLTN data
    ENTLN_file=ENTLN_folder+flash_DATE_STR+'.npy'
    if os.path.exists(ENTLN_file) is True:
        EN=np.load(ENTLN_file)
    else:
        EN=np.array([])
    
    #convert entln time to tod and lat lon to  xy in km
    EN_t_unix_ns=(EN[:,0]*1e9).astype(np.int64)
    EN_tod=((EN_t_unix_ns-ref_t_ns_4_tod)/1e9).reshape(-1,1)
    
    fake_EN_alt=np.zeros(len(EN_tod)).reshape(-1,1) # assume all at zero 
    EN_lat=EN[:,1].reshape(-1,1)
    EN_lon=EN[:,2].reshape(-1,1)
    EN_type=EN[:,3].reshape(-1,1)
    EN_pc=(EN[:,4]/1e3).reshape(-1,1)
    
    EN_lla=np.hstack((EN_lat,EN_lon,fake_EN_alt))
    
    EN_xyz=geodetic_to_enu(EN_lla,ref_lat,ref_lon,ref_alt)
    EN_xy_km=EN_xyz[:,0:2]/1e3
    
    # EN_type are also for fake event ID, just to match format
    EN_events=np.hstack((EN_type,EN_xy_km,fake_EN_alt,EN_tod,EN_pc)) 
    cg_events=EN_events[EN_type.flatten()==0,:]
    
     # keep only flashes with sources more than a threshold
    n_sources_thres=50
    n_sources_col=flash_info[:,2]
        
    # flash info for big flash only
    big_flash_info=flash_info[n_sources_col>n_sources_thres,:]
    big_flash_fidx=big_flash_info[:,3]
    
    #here we keep S_sorted for big flash only
    tf_row_keep=np.isin(S_sorted[:,-1],big_flash_fidx)
    S_sorted_big_flash=S_sorted[tf_row_keep,:]
    
    # big flash info each column: flash_t0, flash_t1, num_sources, flash_idx, flash_type, number of strokes
    big_flash_info,cg_events=flash_type_assign(S_sorted_big_flash,cg_events,big_flash_info)
    
    
    for ii, info in enumerate(big_flash_info):
        print(ii)
        # if ii!=19:
        #     continue
    
        flash_ID=info[3]
        sources_idx=np.where(S_sorted_big_flash[:,-1]==flash_ID)[0]
        
        #check if there is any cg events detected in this flash:
        if info[-1]>0: # if it is a cg flash
            # cg_event format (type,x, y, altitude,tod,en_pc,big_flash_row )
            cg_events_in_this_flash=cg_events[cg_events[:,-1]==ii,:] 
        else:
            cg_events_in_this_flash=np.array([])  
        
        
        # now the f is in the format of (x,y,z,t,chi,power,lla,flash_ID), in which xyz are in km
        f=S_sorted_big_flash[sources_idx]
        
        #lets do a body check to remove noises in lma flash
        sep_thres=1
        f,tf_idx=body_check(f,sep_thres) # tf idx is the boolean mask of kept sources
        
        f_t=f[:,3]
        f_t1=np.min(f_t)
        f_t2=np.max(f_t)

        
        f_t1_hhmmss=tod_2_tstamp(f_t1)[:18]  # we need precision to ns
        f_t1_str=flash_DATE_STR+'T'+f_t1_hhmmss
        f_t1_tstamp=pd.Timestamp(f_t1_str)
        
        f_t2_hhmmss=tod_2_tstamp(f_t2)[:18]  # we need precision to ns
        f_t2_str=flash_DATE_STR+'T'+f_t2_hhmmss
        f_t2_tstamp=pd.Timestamp(f_t2_str)
        
        
        ref_t_stamp=tod_2_tstamp(f_t1)
        ref_t_stamp_2_ms=ref_t_stamp[:12]
        full_t_stamp=flash_DATE_STR+'T'+ref_t_stamp_2_ms
        
        # here we geolocate all pixels at the flash starting time
        all_pixels_coordinates,pxpy=geolocate_all_pixels(t=f_t1_tstamp,one_second_df=l.one_second)
        all_pixels_latlon=np.hstack((all_pixels_coordinates.lat.reshape(-1,1),all_pixels_coordinates.lon.reshape(-1,1)))
        all_pixels_xy=haversine_latlon_xy_conversion(all_pixels_latlon,NALMA_center)
        

        #get fov of LIS at flash t1 and flash t2
        fov1 = ltgLIS.get_fov(l.one_second, times=f_t1_tstamp)
        fov2 = ltgLIS.get_fov(l.one_second, times=f_t2_tstamp)
        
        #we need to make sure that LMA flashes are within the fovs
        fov1_lonlat=np.hstack((fov1[0].lon.reshape(-1,1),fov1[0].lat.reshape(-1,1)))
        fov2_lonlat=np.hstack((fov2[0].lon.reshape(-1,1),fov2[0].lat.reshape(-1,1)))
        
        
        #check if lma sources are all within the fovs
        f_lat=f[:,6]
        f_lon=f[:,7]
        f_latlon=np.hstack((f_lat.reshape(-1,1),f_lon.reshape(-1,1)))
        f_lonlat=np.hstack((f_lon.reshape(-1,1),f_lat.reshape(-1,1)))
        
        f_xy=haversine_latlon_xy_conversion(f_latlon,NALMA_center)
        lma_x=f_xy[:,0]
        lma_y=f_xy[:,1]
        lma_z=f[:,8]/1e3
        lma_t=(f_t-f_t1)*1e3
        
        # make sure all LMA sources are within the fovs of LIS at flash t1 and t2
        in_fov1=check_if_within_polygon(fov1_lonlat,f_lat,f_lon)
        if in_fov1 is False:
            print(ii,'out of FOV1')
            continue 
        
        in_fov2=check_if_within_polygon(fov2_lonlat,f_lat,f_lon)
        if in_fov2 is False:
            print(ii,'out of FOV2')
            continue 
        
        ########################################################################
        ######  if the code get here it means the lma flash is within the FOV of LIS
        ########################################################################
        
        # determine the convex hull of the LMA flash in 2D
        assert len(f_lonlat)>2 # MAKRE SURE lma SOURCES ARE MORE THAN 2
        hull = ConvexHull(f_xy)
        area=hull.volume
        
        ## here we need to get the vertices of the convex hull and based on that create a polygon object
        hull_x=(lma_x[hull.vertices]).reshape(-1,1)
        hull_y=(lma_y[hull.vertices]).reshape(-1,1)
        hull_xy=np.hstack((hull_x,hull_y))
        rr=LinearRing(hull_xy)
        hull_polygon = Polygon(rr)
        hull_polygon_expanded = Polygon(hull_polygon.buffer(2.0).exterior) 
        
        #find the centroid of the LMA polygon (before expanding)
        hull_centroid_x,hull_centroid_y=hull_polygon.centroid.coords.xy
        
        #find the corresponding px and py for the centroid
        d_pixels_2_centroid=d_point_to_points_2d([hull_centroid_x,hull_centroid_y],all_pixels_xy)
        centroid_pixel_idx=np.where((d_pixels_2_centroid<5)&(d_pixels_2_centroid==np.min(d_pixels_2_centroid)))[0][0]
        centroid_pxpy=pxpy[centroid_pixel_idx]
        

        # first lets narrow-down the events data based on temporal and spatial span of lma data        
        #temporal part
        E_tod= group_t_stamps_to_tod(E.time)
        idx_keep_temporal=np.where((E_tod>=f_t1)&(E_tod<=f_t2))[0]
        if len(idx_keep_temporal)==0:
            print(ii,'base on time span of LMA flash, no LIS events found')
            continue
        
        E1=E.iloc[idx_keep_temporal]
        
        #spatial part
        n_lma_sources=len(f_lat)
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
            print(ii,'base on spatial scale of LMA flash, no LIS events found')
            continue
        

        E2=E1.iloc[idx_keep_spatial]
        
        # optional, we could add addtional events that were not in E2 but in parent groups 
        E2=expand_E2(E2,G,E)


        # get total radiance vs time
        t_frame,rad_frame,area_frame= get_radiance_area_vs_time(E2)
        t_frame=(t_frame-f_t1)*1e3
        
        lla = ltgLIS.geolocate_events_in_frame(E2, l.one_second)
        ev_poly = ltgLIS.event_poly(E2, corners=lla, latlon=True, fill=True)
                
        #now we get the ploygon of events in latlon, lets covert latlon of polygon to xy
        ev_poly_latlon,ev_poly_xy = convert_poly_latlon_to_xy(ev_poly)

        fig,axs=LIS_plot_layout()
        
        ax0_xlim=[np.min(lma_t)-10,np.max(lma_t)+10]
        # lma source h vs time plot
        axs[0].scatter(lma_t,lma_z,c=f_t, marker='.', zorder=13,cmap='jet')
        #plot cg events if any:
        if len(cg_events_in_this_flash)>0:
            ax0_ylim=axs[0].get_ylim()
            for cg in cg_events_in_this_flash:
                cg_t=(cg[4]-f_t1)*1e3
                axs[0].plot([cg_t,cg_t],[ax0_ylim[0],ax0_ylim[1]],'-k')
            axs[0].set_ylim(ax0_ylim[0],ax0_ylim[1])
        axs[0].set_xlim(ax0_xlim)
        
        # radiance/area vs time plot
        axs[1].plot(t_frame,rad_frame,'r.')
        axs[1].bar(t_frame,rad_frame,2,facecolor='r')
        axs[1].set_xlim(ax0_xlim)
        axs[1].set_xlabel('Milliseconds after '+ full_t_stamp)
        
        E2_radiance=E2.radiance.values
        for m, poly in enumerate(ev_poly_xy):
            # p = Polygon(poly, facecolor = 'k')
            # ax.add_patch(p)
            rad=E2_radiance[m]
            
            if rad<3500:
                rad=3500
            if rad>40000:
                rad=40000
            
            color_scale=(rad-3500)/(40000-3500) # this value range from 0 to 1
            axs[2].fill(poly[:,0],poly[:,1],facecolor=cmap(color_scale), edgecolor='black', linewidth=1)
        axs[2].scatter(lma_x,lma_y,c=f_t, marker='.', zorder=13,cmap='jet')
                
        # get the xy bounds of the expanded polygon, needs this for deteming plot limit
        ep_minx, ep_miny, ep_maxx, ep_maxy=hull_polygon_expanded.bounds
        
        
        # plot LMA flash centroid
        axs[2].scatter(hull_centroid_x,hull_centroid_y, c='k', s=150, marker='X', zorder=14)
                
        # Now we can check if the expanded hull_polygon intersect with any pixels in E2
        LIS_LMA_intesect_tf= LIS_LMA_intesect(ev_poly_xy,hull_polygon_expanded)
        if LIS_LMA_intesect_tf ==True:
            print('DETECTED')
        else:
            print('MISSED')

                
        x1, y1 = hull_polygon.exterior.xy
        x2, y2 = hull_polygon_expanded.exterior.xy
        
        axs[2].plot(x1,y1)
        axs[2].plot(x2,y2)
        
        # here we need to make the plot square, to do this, we need to make sure
        ev_poly_xy_ary=np.concatenate(ev_poly_xy)
        all_x=np.concatenate((ev_poly_xy_ary[:,0],lma_x,np.array([ep_minx, ep_maxx])))
        all_y=np.concatenate((ev_poly_xy_ary[:,1],lma_y,np.array([ep_miny, ep_maxy])))
        
        min_x=np.min(all_x)
        max_x=np.max(all_x)
        min_y=np.min(all_y)
        max_y=np.max(all_y)
        
        #length compensation
        length_x=max_x-min_x
        length_y=max_y-min_y
        
        if length_x>length_y:
            extra=(length_x-length_y)/2
            min_y=min_y-extra
            max_y=max_y+extra
        
        if length_x<length_y:
            extra=(length_y-length_x)/2
            min_x=min_x-extra
            max_x=max_x+extra
        
        
        pad_x=(max_x-min_x)*0.1
        pad_y=(max_y-min_y)*0.1
        
        ax2_xlim=[min_x-pad_x,max_x+pad_x]
        ax2_ylim=[min_y-pad_y,max_y+pad_y]
        
        
        axs[2].set_xlim(ax2_xlim)
        axs[2].set_ylim(ax2_ylim)
                    
        axs[2].set_aspect('equal', adjustable='box')
        plt.show()
        
         
        # get rid of puncuation that can not be used as 
        full_t_stamp=full_t_stamp.replace(':','_')
        full_t_stamp=full_t_stamp.replace('.','_')
        fig_name=fig_folder+str(ii)+'_'+full_t_stamp+'.png'
        
        fig.savefig(fig_name,dpi=300,bbox_inches='tight')

      
