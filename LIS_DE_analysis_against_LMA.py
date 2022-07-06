import numpy as np
import matplotlib.pyplot as plt
import pickle

import pandas as pd

from matplotlib import cm
from matplotlib.gridspec import GridSpec
from scipy.spatial import ConvexHull
from shapely.geometry import LinearRing
from shapely.geometry.polygon import Polygon

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

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

def get_limits_for_same_scale_planview (hull_polygon_expanded,lma_x,lma_y,ev_poly_xy=[]):
    
    
    ep_minx, ep_miny, ep_maxx, ep_maxy=hull_polygon_expanded.bounds
    
    if len(ev_poly_xy)!=0:
        ev_poly_xy_ary=np.concatenate(ev_poly_xy)
        all_x=np.concatenate((ev_poly_xy_ary[:,0],lma_x,np.array([ep_minx, ep_maxx])))
        all_y=np.concatenate((ev_poly_xy_ary[:,1],lma_y,np.array([ep_miny, ep_maxy])))
    else:
        all_x=np.concatenate((lma_x,np.array([ep_minx, ep_maxx])))
        all_y=np.concatenate((lma_y,np.array([ep_miny, ep_maxy])))
    
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
    
    return ax2_xlim, ax2_ylim    

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

def plot_LMA_LIS_MATCH(m,ii,LMA_center,fig_save=False,fig_folder=None):
    
    #set the colormap
    cmap=cm.get_cmap('binary',100)
    
    
    # LIS-LMA match plot
    fig,axs=LIS_LMA_plot_layout()
    
    f=m['LMA']

    f_t=(f['t'].astype('int64')/int(1e9)).to_numpy()
    f_t1=np.min(f_t)
    
    f_lat=f['lat'].to_numpy()
    f_lon=f['lon'].to_numpy()
    f_latlon=np.hstack((f_lat.reshape(-1,1),f_lon.reshape(-1,1)))
    f_xy=haversine_latlon_xy_conversion(f_latlon,LMA_center)
    lma_x=f_xy[:,0]
    lma_y=f_xy[:,1]
    
    f_t1_stamp=str(pd.to_datetime(f_t1, unit='s', origin='unix'))
    
    lma_t=(f_t-f_t1)*1e3 # convert it to ms
    ax0_xlim=[np.min(lma_t)-10,np.max(lma_t)+10]
    
    lma_z=f['altitude'].to_numpy()
    
    # get the xy bounds of the expanded polygon, needs this for deteming plot limit
    # get the convexhull polygon of the lma flash
    # TODO: here need a line for hull in lat and lon
    hull_polygon=get_lma_convex_hull_polygon(f_xy)
    
    #find the centroid of the LMA polygon (before expanding)
    hull_centroid_x,hull_centroid_y=hull_polygon.centroid.coords.xy
    
    # expand the convex hull by 2 km to account of location offset of LIS
    hull_polygon_expanded = Polygon(hull_polygon.buffer(2.0).exterior) 
    
    x1, y1 = hull_polygon.exterior.xy
    x2, y2 = hull_polygon_expanded.exterior.xy
    
    ######################################
    # lma source height vs time plot
    ######################################
    axs[0].scatter(lma_t,lma_z,c=f_t, marker='.', zorder=13,cmap='jet')
    
    #plot cg events if any:
    if 'RS' in m.keys():
        ax0_ylim=axs[0].get_ylim()
        cg_events_in_this_flash=m['RS']
        cg_events_in_this_flash['t']=cg_events_in_this_flash['t'].astype('int64')/1e9
        cg_events_in_this_flash=cg_events_in_this_flash.to_numpy()
        
        for cg in cg_events_in_this_flash:
            cg_t=(cg[3]-f_t1)*1e3
            axs[0].plot([cg_t,cg_t],[ax0_ylim[0],ax0_ylim[1]],'-k')
        axs[0].set_ylim(ax0_ylim[0],ax0_ylim[1])
    axs[0].set_xlim(ax0_xlim)
    
    
    ######################################
    # LIS total radiance vs time plot
    #####################################
    if m['LIS detection']==True:
        # get total radiance vs time
        E2=m['LIS_events']
        t_frame,rad_frame,area_frame= get_radiance_area_vs_time(E2)
        t_frame=(t_frame-f_t1)*1e3
              
        ev_poly_xy = m['LIS_events_polygon_xy']
        
        # radiance/area vs time plot
        axs[1].plot(t_frame,rad_frame,'r.')
        axs[1].bar(t_frame,rad_frame,2,facecolor='r')
        axs[1].set_xlim(ax0_xlim)
        axs[1].set_xlabel('Milliseconds after '+ f_t1_stamp)
        
        ###############################################
        # plan view  lma overlaid on LIS events polygon
        ###############################################
        
        E2_radiance=E2.radiance.values
        for jj, poly in enumerate(ev_poly_xy):
            # p = Polygon(poly, facecolor = 'k')
            # ax.add_patch(p)
            rad=E2_radiance[jj]
            
            if rad<3500:
                rad=3500
            if rad>40000:
                rad=40000
            
            color_scale=(rad-3500)/(40000-3500) # this value range from 0 to 1
            axs[2].fill(poly[:,0],poly[:,1],facecolor=cmap(color_scale), edgecolor='black', alpha=0.5, linewidth=2)
        
    # lma plan view plot
    axs[2].scatter(lma_x,lma_y,c=f_t, marker='.', zorder=13,cmap='jet')
    # plot LMA flash centroid
    axs[2].scatter(hull_centroid_x,hull_centroid_y, c='k', s=150, marker='X', zorder=14)
    # hull polygons, original and expanded
    axs[2].plot(x1,y1)
    axs[2].plot(x2,y2)
    
    if m['LIS detection']==True:
        ax2_xlim, ax2_ylim = get_limits_for_same_scale_planview (hull_polygon_expanded,lma_x,lma_y,ev_poly_xy)
    else:
        ax2_xlim, ax2_ylim = get_limits_for_same_scale_planview (hull_polygon_expanded,lma_x,lma_y)
        
    
    axs[2].set_xlim(ax2_xlim)
    axs[2].set_ylim(ax2_ylim)
    
    # get rid of puncuation that can not be used as 
    if fig_save:
        full_t_stamp=f_t1_stamp.replace(':','_')
        full_t_stamp=full_t_stamp.replace('.','_')
        fig_name=fig_folder+str(ii)+'_'+full_t_stamp+'.png'
        fig.savefig(fig_name,dpi=300,bbox_inches='tight')



M=load_obj('E:/NALMA_LIS/NALMA_LIS_matches.pkl')
LMA_center=np.array([34.8,-86.85])


# 0:DE 
# 1:Flash area
# 2:Number of sources
# 3:Median source altitude
# 4:Flash type
# 5:max EN CG peak current
# 6：Day or night
# 7:px
# 8:py

F=np.zeros((len(M),9))


for i, m in M.items():
        
    f=m['LMA']
    
    # if LIS detected 
    if m['LIS detection'] == True:
        F[i,0]=1
    else:
        F[i,0]=0
    
    # flash area
    F[i,1]=np.ceil(m['flash area'])
    
    # Number of sources
    F[i,2]=len(f)
    
    # Median source height
    F[i,3]=np.median(f['altitude'])
    
    # Type
    if m['type']=='IC flash':
        F[i,4]=0
    elif m['type']=='Negative CG flash':
        F[i,4]=1
    elif m['type']=='Positive CG flash':
        F[i,4]=2
    else:
        F[i,4]=3 # bipolar flash
        
    # maximum CG peak current in absoluate value
    if 'RS' in m.keys():
        F[i,5]=np.max(abs(m['RS']['peak current']))
    else:
        F[i,5]=-999
        
    # day or night
    if m['dn']=='day':
        F[i,6]=1
    else:
        F[i,6]=0
        
    # px
    F[i,7]=m['centroid pxpy'][0]
        
    # py
    F[i,8]=m['centroid pxpy'][1]
    

# Calculte DE
num_LMA_flashes=len(F)
num_LIS_detections=int(np.sum(F[:,0]))
DE=np.around(num_LIS_detections/num_LMA_flashes,2)
print(f"LIS detected {num_LIS_detections}/{num_LMA_flashes} LMA flashes, DE is {DE}")

# Calculte day/night DE
F_day=F[F[:,6]==1,:]
F_night=F[F[:,6]==0,:]
num_LMA_flashes_day=len(F_day)
num_LMA_flashes_night=len(F_night)

num_LIS_detections_day=int(np.sum(F_day[:,0]))
num_LIS_detections_night=int(np.sum(F_night[:,0]))

DE_day=np.around(num_LIS_detections_day/num_LMA_flashes_day,2)
DE_night=np.around(num_LIS_detections_night/num_LMA_flashes_night,2)

print(f"During day: LIS detected {num_LIS_detections_day}/{num_LMA_flashes_day} LMA flashes, DE is {DE_day}")
print(f"During day: LIS detected {num_LIS_detections_night}/{num_LMA_flashes_night} LMA flashes, DE is {DE_night}")

#########################################
#plot DE with different variables
#########################################

DE_col=F[:,0]
FA_col=F[:,1]
SA_col=F[:,3]

# flash_area histogram
fig,ax=plt.subplots(figsize=(10,8))
bins=2**(np.arange(13)) # TODO: NEED TO CALCULATE THE ORDER

bins_middle=np.diff(bins)/2+bins[0:-1]

ax.hist(FA_col,bins=bins)
ax.set_xscale('log')
ax.set_xlabel('Flash area ($km^2$)')
ax.set_ylabel('Number')

# DE vs flash area
bin_idx=np.digitize(F[:,1],bins=bins)
DE_per_bin=np.zeros(len(bins)-1)

bin_idx_unique=np.arange(len(bins)-1)+1
for i, b_idx in enumerate(bin_idx_unique):
    idx=np.where(bin_idx==b_idx)[0]
    DE_bin=DE_col[idx]
    DE_percent_bin=np.sum(DE_bin)/len(DE_bin)
    DE_per_bin[i]=DE_percent_bin
    
ax1=ax.twinx()
ax1.plot(bins_middle,DE_per_bin,'-X',color='r')
ax1.set_ylabel('Detection Efficiency')


# Meidan source altitude histogram
fig,ax=plt.subplots(figsize=(10,8))
bins=np.arange(3,18,0.5) 
# remove flashes with altititude that are likely wrong
idx_keep=np.where((SA_col>=3)&(SA_col<=18))[0]
SA_col=SA_col[idx_keep]

bins_middle=np.diff(bins)/2+bins[0:-1]

ax.hist(SA_col,bins=bins)
ax.set_xlabel('Median Source Altitude (km)')
ax.set_ylabel('Number')

# DE vs source altitude
bin_idx=np.digitize(SA_col,bins=bins)
DE_per_bin=np.zeros(len(bins)-1)

bin_idx_unique=np.arange(len(bins)-1)+1
for i, b_idx in enumerate(bin_idx_unique):
    idx=np.where(bin_idx==b_idx)[0]
    DE_bin=DE_col[idx]
    DE_percent_bin=np.sum(DE_bin)/len(DE_bin)
    DE_per_bin[i]=DE_percent_bin
    
ax1=ax.twinx()
ax1.plot(bins_middle,DE_per_bin,'-X',color='r')
ax1.set_ylabel('Detection Efficiency')


# num_LIS_detection=num_LIS_detection_day+num_LIS_detection_night

# DE=np.around(num_LIS_detection/num_flashes,2)
# DE_day=np.around(num_LIS_detection_day/num_day_flashes,2)
# DE_night=np.around(num_LIS_detection_night/num_night_flashes,2)


# print(f"LIS detected {num_LIS_detection}/{len(M)} LMA flashes, DE is {DE}")
# print(f"During day: LIS detected {num_LIS_detection_day}/{num_day_flashes} LMA flashes, DE is {DE_day}")
# print(f"During night: LIS detected {num_LIS_detection_night}/{num_night_flashes} LMA flashes, DE is {DE_night}")



# ii=29
# m=M[ii]
# plot_LMA_LIS_MATCH(m,ii,LMA_center)