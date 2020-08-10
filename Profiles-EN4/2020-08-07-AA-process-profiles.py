#!/usr/bin/env python
# coding: utf-8

# ## Import all the libraries

# In[ ]:


import numpy as np
import dask
import xarray as xr
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import glob as glob
import time
from dask.diagnostics import ProgressBar
from datetime import date
import json
import os
import warnings
import seawater
import re

warnings.filterwarnings('ignore')


# ## Parameters

# In[ ]:


ymin=2010;mmin=1;dmin=1
ymax=2010;mmax=9;dmax=30
# depth of the desired comparison profile in m
depthmin=1000
# radius of the circle around the profile location in which we take the modeled profiles, in °  
radius_max=0.25
# period of time around the profile sampling date in which we take the modeled profiles, in days
period=5
# minimum amount of model profiles to be considered to make a significant statistical comparison, for instance in a 1° square and 30-days window we have 2.6 millions modeled profiles, in a 0.5°x10 days 216 000
number_of_model_profiles=100000
# 
debug_plot=False


# In[ ]:


datemin=datetime.date(ymin,mmin,dmin)
datemax=datetime.date(ymax,mmax,dmax)
jsonfile='txt/MEDWEST60-BLBT02_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'


# ## Read the jsonfile

# In[ ]:


sourcefile=open(jsonfile,'rU')
infos=json.load(sourcefile)
nb_profilesEN4=len(infos)


# ## Loop on the number of profiles

# In[ ]:


def process_one_profile(prof):
    # Process one profile
    list_profiles = infos.keys()
    reference =  str(list(list_profiles)[prof])
    print('Processing profile ', reference)

    # Get all infos from json file
    lat_prof = infos[list(list_profiles)[prof]]['latitude']
    lon_prof = infos[list(list_profiles)[prof]]['longitude']
    date_prof = infos[list(list_profiles)[prof]]['date']
    file_prof = infos[list(list_profiles)[prof]]['file']
    prof_prof = infos[list(list_profiles)[prof]]['profile no']
    i0 = infos[list(list_profiles)[prof]]['i0']
    j0 = infos[list(list_profiles)[prof]]['j0']

    # List of all model files involved
    date_profmin=datetime.date(int(date_prof[0:4]),int(date_prof[5:7]),int(date_prof[8:10]))-datetime.timedelta(days=int(period))
    date_profmax=datetime.date(int(date_prof[0:4]),int(date_prof[5:7]),int(date_prof[8:10]))+datetime.timedelta(days=int(period))
    def date_range(start, end):
        r = (end+datetime.timedelta(days=1)-start).days
        return [start+datetime.timedelta(days=i) for i in range(r)]
    dateList = date_range(date_profmin, date_profmax) 
    dirmod="/gpfsstore/rech/egi/commun/MEDWEST60/extracted_eNATL60/allv/"
    list_filesmod_T=[]
    list_filesmod_S=[]
    for date in dateList:
        year=date.year
        month=date.month
        day=date.day
        mm="{:02d}".format(month) #month on 2 digits
        dd="{:02d}".format(day) # day on 2 digits
        list_filesmod_T.append(dirmod+'MEDWEST60-BLBT02_y'+str(year)+'m'+str(mm)+'d'+str(dd)+'.1h_gridT.nc')
        list_filesmod_S.append(dirmod+'MEDWEST60-BLBT02_y'+str(year)+'m'+str(mm)+'d'+str(dd)+'.1h_gridS.nc')
    print(list_filesmod_T)  
    
    # Read model files
    dsT=xr.open_mfdataset(list_filesmod_T)
    dsS=xr.open_mfdataset(list_filesmod_S)
    gdpts=np.int(np.round(radius_max*60))

    lonmod=dsT.nav_lon[0,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    latmod=dsT.nav_lat[0,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    tempmod=dsT.votemper[:,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    saltmod=dsS.vosaline[:,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    depthmod=dsT.deptht
    meshfile='/gpfsstore/rech/egi/commun/MEDWEST60/MEDWEST60-I/mesh_mask.nc'
    ds=xr.open_dataset(meshfile)
    maskmod=ds.tmask[0,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    maskmod0=ds.tmask[0,0,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]

    #Get the depth at every grid point
    d,ly,lx=maskmod.shape
    depthmod2d=np.zeros([lx,ly])
    for j in np.arange(ly):
        for i in np.arange(lx):
            depthmod2d[j,i]=depthmod[np.min(np.where(maskmod[:,j,i].values<1))].values
    
    # Read the profile file
    diren4="/gpfswork/rech/egi/rote001/EN4/"
    tfileEN4=diren4+file_prof
    dsen4=xr.open_dataset(tfileEN4)
    laten4=dsen4['LATITUDE'][prof_prof]
    lonen4=dsen4['LONGITUDE'][prof_prof]
    dayen4=dsen4['JULD'][prof_prof]
    tempen4=dsen4['TEMP'][prof_prof]
    salten4=dsen4['PSAL'][prof_prof]
    presen4=dsen4['PRES_ADJUSTED'][prof_prof]
    depthen4=seawater.dpth(presen4,laten4)

    # Find the last level and reduce the profiles
    indzprof=np.min(np.where(np.isnan(depthen4)==True))
    dmax=depthen4[indzprof-1]
    obsred_dep=np.zeros(int(indzprof))
    obsred_temp=np.zeros(int(indzprof))
    obsred_salt=np.zeros(int(indzprof))
    for z in np.arange(int(indzprof)):
        obsred_dep[int(z)]=depthen4[int(z)]
        obsred_temp[int(z)]=tempen4[int(z)]
        obsred_salt[int(z)]=salten4[int(z)]

    # Find the model profiles
    lon_stacked = lonmod.stack(profile=('x', 'y'))
    lat_stacked = latmod.stack(profile=('x', 'y'))
    mask_stacked = maskmod0.stack(profile=('x', 'y'))
    xr_depthmod2d=xr.DataArray(depthmod2d, dims=("x", "y"))
    depth_stacked = xr_depthmod2d.stack(profile=('x', 'y'))
    
    distance_threshold = radius_max
    square_distance_to_observation = (lon_stacked - lon_prof)**2 + (lat_stacked-lat_prof)**2
    square_distance_to_observation_mask = np.ma.masked_where(mask_stacked==0.,square_distance_to_observation) 
    square_distance_to_observation_sorted = np.sort(square_distance_to_observation_mask)
    nb_profiles_per_timestep=number_of_model_profiles/(24*period*2+24)
    new_threshold=square_distance_to_observation_sorted[int(np.round(nb_profiles_per_timestep)+1)]
    is_closer_to_observation = (square_distance_to_observation < new_threshold) & (depth_stacked > depthmin)

    model_temperature_stacked = tempmod.stack(profile=('x', 'y'))
    model_salinity_stacked = saltmod.stack(profile=('x', 'y'))

    model_temperature_near_observation = model_temperature_stacked.where(is_closer_to_observation,drop=True)
    model_salinity_near_observation = model_salinity_stacked.where(is_closer_to_observation, drop=True)
    lat_near_observation = lat_stacked.where(is_closer_to_observation, drop=True)
    lon_near_observation = lon_stacked.where(is_closer_to_observation, drop=True)
    
    # Compute statistics on the model profiles
    temp_model_mean = np.mean(model_temperature_near_observation,axis=(0,2))
    temp_percentile_10= np.percentile(model_temperature_near_observation,10,axis=(0,2))
    temp_percentile_90= np.percentile(model_temperature_near_observation,90,axis=(0,2))
    salt_model_mean = np.mean(model_salinity_near_observation,axis=(0,2))
    salt_percentile_10= np.percentile(model_salinity_near_observation,10,axis=(0,2))
    salt_percentile_90= np.percentile(model_salinity_near_observation,90,axis=(0,2))
    
    # Interpolate on obs vertical grid
    temp_model_mean_depobs=np.interp(obsred_dep,depthmod,temp_model_mean)
    temp_model_percentile_10_depobs=np.interp(obsred_dep,depthmod,temp_percentile_10)
    temp_model_percentile_90_depobs=np.interp(obsred_dep,depthmod,temp_percentile_90)
    salt_model_mean_depobs=np.interp(obsred_dep,depthmod,salt_model_mean)
    salt_model_percentile_10_depobs=np.interp(obsred_dep,depthmod,salt_percentile_10)
    salt_model_percentile_90_depobs=np.interp(obsred_dep,depthmod,salt_percentile_90)
    
    # Make a debug plot
    if debug_plot == True:
        fig, axs = plt.subplots(1,2, figsize=(10, 6))
        axs = axs.ravel()
        title = 'Temperature and Salinity Profiles for profile '+reference
        plt.suptitle(title,size = 25,y=1.05)
        axs[0].plot(temp_model_mean_depobs,obsred_dep,'b.-', label='temp model')
        axs[0].plot(obsred_temp,obsred_dep,'k.-', label='temp en4')
        axs[0].set_ylabel('Depth [m]', size=14)
        axs[0].set_ylim(2000, 0)
        axs[0].grid(True, which='both')
        axs[0].xaxis.tick_top()
        axs[0].xaxis.set_label_position('top') 
        axs[0].plot(temp_model_percentile_10_depobs,obsred_dep,'b-', label='percent10')
        axs[0].plot(temp_model_percentile_90_depobs,obsred_dep,'b-', label='percent90')
        axs[0].fill_betweenx(obsred_dep, temp_model_percentile_10_depobs, x2=temp_model_percentile_90_depobs, alpha=0.2, facecolor='b')

        axs[1].plot(salt_model_mean_depobs,obsred_dep,'b.-', label='salt model')
        axs[1].plot(obsred_salt,obsred_dep,'k.-', label='salt en4')
        axs[1].set_ylabel('Depth [m]', size=14)
        axs[1].set_ylim(2000, 0)
        axs[1].grid(True, which='both')
        axs[1].xaxis.tick_top()
        axs[1].xaxis.set_label_position('top') 
        axs[1].plot(salt_model_percentile_10_depobs,obsred_dep,'b-', label='percent10')
        axs[1].plot(salt_model_percentile_90_depobs,obsred_dep,'b-', label='percent90')
        axs[1].fill_betweenx(obsred_dep, salt_model_percentile_10_depobs, x2=salt_model_percentile_90_depobs, alpha=0.2, facecolor='b')
        fig.tight_layout()
        plt.savefig('figs/MEDWEST60-BLBT02_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'_prof'+str(prof)+'.png')

    # Write netcdf file
    match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
    debut_ref=match.group(1)
    fin_ref=match.group(2)
    dirname=diren4+'profiles_files/MEDWEST60-BLBT02/'
    if not os.path.exists(dirname):
        os.makedirs(dirname)    
  
    outname=dirname+str(debut_ref)+str(fin_ref)+'_MEDWEST60-BLBT02_TS.nc'
    print('output file is '+outname)
    dsout=Dataset(outname,'w')

    today=date.today()
    dsout.description = "This file contains one profile of temperature and salinity from EN4 dataset and the mean and 10 and 90 percentile of NATL60-CJM165 data within a 0.25deg circle around the location of the profile and 15 days before and after it has been sampled. This file has been created "+str(today.day)+"/"+str(today.month)+"/"+str(today.year)

    depth=dsout.createDimension('depth',len(obsred_dep))
    x=dsout.createDimension('x',1)
    y=dsout.createDimension('y',1)
    
    lat = dsout.createVariable('latitude_profileEN4', 'f8', ('y','x'))
    lat.standart_name="latitude_profileEN4"
    lat.long_name = "Latitude of selected EN4 profile" 
    lat.units = "degrees_north"

    lon = dsout.createVariable('longitude_profileEN4', 'f8', ('y','x'))
    lon.standart_name="longitude_profileEN4"
    lon.long_name = "Longitude of selected EN4 profile" 
    lon.units = "degrees_east"

    time = dsout.createVariable('time_profileEN4', 'f8', ('y','x'))
    time.standart_name="time_profileEN4"
    time.timeg_name = "Time in seconds from 1-1-1958 of selected EN4 profile" 
    time.units = "s"

    depth_en4 = dsout.createVariable('depth_en4', 'f8', ('depth'),fill_value=0.)
    depth_en4.units = "m" 
    depth_en4.valid_min = 0.
    depth_en4.valid_max = 8000.
    depth_en4.long_name = "Depth" 

    temp_en4 = dsout.createVariable('temp_profileEN4', 'f8', ('depth'),fill_value=0.)
    temp_en4.units = "degC" 
    temp_en4.valid_min = -10.
    temp_en4.valid_max = 40.
    temp_en4.long_name = "Temperature profile of the selected EN4 profile" 

    salt_en4 = dsout.createVariable('salt_profileEN4', 'f8', ('depth'),fill_value=0.)
    salt_en4.units = "PSU" 
    salt_en4.valid_min = 20.
    salt_en4.valid_max = 40.
    salt_en4.long_name = "Salinity profile of the selected EN4 profile" 

    mean_temp_model = dsout.createVariable('mean_temp_model', 'f8', ('depth'),fill_value=0.)
    mean_temp_model.units = "degC" 
    mean_temp_model.valid_min = -10.
    mean_temp_model.valid_max = 40.
    mean_temp_model.long_name = "Mean Temperature profile of the model" 

    mean_salt_model = dsout.createVariable('mean_salt_model', 'f8', ('depth'),fill_value=0.)
    mean_salt_model.units = "PSU" 
    mean_salt_model.valid_min = 20.
    mean_salt_model.valid_max = 40.
    mean_salt_model.long_name = "Mean Salinity profile of the model" 

    percent10_temp_model = dsout.createVariable('percent10_temp_model', 'f8', ('depth'),fill_value=0.)
    percent10_temp_model.units = "degC" 
    percent10_temp_model.valid_min = -10.
    percent10_temp_model.valid_max = 40.
    percent10_temp_model.long_name = "Percent 10 Temperature profile of the model" 

    percent10_salt_model = dsout.createVariable('percent10_salt_model', 'f8', ('depth'),fill_value=0.)
    percent10_salt_model.units = "PSU" 
    percent10_salt_model.valid_min = 20.
    percent10_salt_model.valid_max = 40.
    percent10_salt_model.long_name = "Percent 10 Salinity profile of the model" 

    percent90_temp_model = dsout.createVariable('percent90_temp_model', 'f8', ('depth'),fill_value=0.)
    percent90_temp_model.units = "degC" 
    percent90_temp_model.valid_min = -90.
    percent90_temp_model.valid_max = 40.
    percent90_temp_model.long_name = "Percent 90 Temperature profile of the model" 

    percent90_salt_model = dsout.createVariable('percent90_salt_model', 'f8', ('depth'),fill_value=0.)
    percent90_salt_model.units = "PSU" 
    percent90_salt_model.valid_min = 20.
    percent90_salt_model.valid_max = 40.
    percent90_salt_model.long_name = "Percent 90 Salinity profile of the model" 


    lat[:]=lat_prof
    lon[:]=lon_prof
    time[:]=(datetime.datetime(int(date_prof[0:4]),int(date_prof[5:7]),int(date_prof[8:10]))-datetime.datetime(1958,1,1,0,0)).total_seconds()
    depth_en4[:]=obsred_dep
    temp_en4[:]=obsred_temp
    salt_en4[:]=obsred_salt
    mean_temp_model[:]=temp_model_mean_depobs
    mean_salt_model[:]=salt_model_mean_depobs
    percent10_temp_model[:]=temp_model_percentile_10_depobs
    percent10_salt_model[:]=salt_model_percentile_10_depobs
    percent90_temp_model[:]=temp_model_percentile_90_depobs
    percent90_salt_model[:]=salt_model_percentile_90_depobs
    dsout.close()  # close the new file
    


# In[ ]:


print("Nb de profiles : "+str(nb_profilesEN4))
print(time.strftime('%d/%m/%y %H:%M',time.localtime()))

for prof in range(nb_profilesEN4):
    list_profiles = infos.keys()
    reference =  str(list(list_profiles)[prof])
    match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
    debut_ref=match.group(1)
    fin_ref=match.group(2)
    diren4="/gpfswork/rech/egi/rote001/EN4/"
    dirname=diren4+'profiles_files/MEDWEST60-BLBT02/'
    outname=dirname+str(debut_ref)+str(fin_ref)+'_MEDWEST60-BLBT02_TS.nc'
    if not os.path.exists(outname):
        process_one_profile(prof)

print(time.strftime('%d/%m/%y %H:%M',time.localtime()))


# In[ ]:




