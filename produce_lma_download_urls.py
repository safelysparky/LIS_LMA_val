# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 08:48:50 2022

This script serves to produce a downloading urls of LMA files, based on LIS_LMA_matched_filenames.txt files.
The ouput url txt file can be used for wget command to download lma files so that you do not need to manually
find corresponding lma files and then download them.

But you need to provide the example url links with placeholders, sometimes LMA filename format 
varies, so some adjustment are needed to match the format


@author: yanan
"""

import pandas as pd
import numpy as np


strftime_format1='%y%m%d_%H%M%S' # e.g. '220412_144000'
strftime_format2='%Y%m%d_%H%M%S' # e.g. '20220412_144000'
strftime_format2='%Y%m%d_%H' # e.g. '20220412_14'



LMA_name='NALMA'
cvs_fname='C:/Users/yanan/Desktop/git_local/LIS_LMA_val/LIS_NALMA_matched_filenames.txt'

one_lma_file_duration=600

one_lma_file_duration_str=str(one_lma_file_duration)
#pad a zero before to make 4 digits
if len(one_lma_file_duration_str)<=3:
    one_lma_file_duration_str='0'+one_lma_file_duration_str

#url_example1="https://data.ghrc.earthdata.nasa.gov/ghrcw-protected/nalma__1/NALMA_220619_234000_0600.dat.gz"
#url_example2="https://data.ghrc.earthdata.nasa.gov/ghrcw-protected/goesrpltnalma__1/goesr_plt_NALMA_20170531_23_3600.dat.gz"
#url_example3="https://data.ghrc.earthdata.nasa.gov/ghrcw-protected/goesrpltksclma__1/goesr_plt_KSCLMA_20170531_235000_0600.dat.gz"
#url_example4="https://data.ghrc.earthdata.nasa.gov/ghrcw-protected/goesrpltcolma__1/goesr_plt_COLMA_20170301_000000_0600.dat.gz"
#url_example5="https://data.nssl.noaa.gov/thredds/fileServer/WRDD/OKLMA/2018/02/24/LYLOUT_180224_000000_0600.dat.gz"



url_with_placeholders="https://data.ghrc.earthdata.nasa.gov/ghrcw-protected/nalma__1/NALMA_{YYMMDD_mmhhss}_{duration}.dat.gz"

strftime_format=strftime_format1

df=pd.read_csv(cvs_fname,names=['LIS_fnames','num_of_flashes','t1','t2'])


LMA_files_url=open(LMA_name+"_files_url.txt", "w")
for index, row in df.iterrows():
    t1_close_flashes=pd.Timestamp(row['t1'])-pd.Timedelta(3, "seconds") # add 3 second buffer
    t2_close_flashes=pd.Timestamp(row['t2'])+pd.Timedelta(3, "seconds") # add 3 second buffer
    
    
    
    t1_epoch_s=np.floor(t1_close_flashes.value/int(one_lma_file_duration*1e9))*one_lma_file_duration
    t2_epoch_s=np.floor(t2_close_flashes.value/int(one_lma_file_duration*1e9))*one_lma_file_duration
        
    ts_epoch_s=np.arange(t1_epoch_s,t2_epoch_s+one_lma_file_duration,one_lma_file_duration)
    
    
    for t_epoch in ts_epoch_s:
        t_epoch_stamp=pd.to_datetime(t_epoch, unit='s', origin='unix')
        # example: "170301_191000" strptime to some format, might need changes depends on the names of lma files
        date_time_str=t_epoch_stamp.strftime(strftime_format) 
        url=url_with_placeholders.format(YYMMDD_mmhhss=date_time_str,duration=one_lma_file_duration_str)# fill the placeholders
        LMA_files_url.write(f"{url}\n")
        
LMA_files_url.close()