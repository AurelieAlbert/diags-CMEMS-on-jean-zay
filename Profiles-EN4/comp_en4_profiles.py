# comp_en4_profiles.py
'''
 Selection of EN4 profiles to compare to a simulation
 In this script the EN4 database is browsed in order to find the profiles that will be compare to the outputs of one simulation. 
 The comparison will not be a simple colocation of the profile inside the model grid but we want to make a statistical comparison of the observed profile with a significant number of profiles close to it in the model (close in terms of space and time, for instance in a 0.5° radius around the profile location and 10 days before and after it has been sampled)
 Therefore the selected profiles must have to be relevant according to some criteria :
   - they must be inside the domain of simulation (the mask file of the configuration file must be provided)  
   - they must be sampled inside the period of simulation (the period shortened by a certain amount of days must be provided)
   - they must go as deep as a given depth (according to the desired depth for the comparison profiles)
   - they must be in a location where there are enough model profiles (for instance if profile is too close to an island or the coast)
'''
# Imports of librairies
import numpy as np

import time
import sys
import os
import glob
import io

import datetime
from datetime import date

from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import dask

import json
import re

import seawater
import cmocean

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import warnings
warnings.filterwarnings('ignore')

# Define the criteria for one profile

# Localisation of the profile location inside model domain
def loc(fileEN4,prof):
    # open the maskfile and get lat lon and mask
    ds=xr.open_dataset(meshfile)
    lat=ds.nav_lat
    lon=ds.nav_lon
    latmin,latmax,lonmin,lonmax=(lat.min(),lat.max(),lon.min(),lon.max())
    tmask=ds.tmask    
    #open the profile file and read the infos on latitude, longitude and date 
    tfileEN4=diren4+fileEN4
    dsen4=xr.open_dataset(tfileEN4)
    laten4=dsen4['LATITUDE'][prof].values
    lonen4=dsen4['LONGITUDE'][prof].values
    
    if (lonmin.values<lonen4)&(lonen4<lonmax.values)&(latmin.values<laten4)&(laten4<latmax.values):
        if not os.path.exists('txt'):
            os.makedirs('txt')    
        with open('txt/prof0.txt','w') as txt_file:
            txt_file.write('Profile_'+str(fileEN4)+'_'+str(prof)+' '+str(lonen4)+' '+str(laten4))
            if os.path.exists('ij_found.out'):
                get_ipython().system('rm ij_found.out')
            
            os.system('/gpfswork/rech/egi/rote001/git/sosie/bin/ij_from_lon_lat.x -i $meshfile -p txt/prof0.txt > txt/output')

            with open('ij_found.out','r') as txt_file:
                last_line = txt_file.readlines()[-1]
                if last_line[1] == '#':
                    print('Profile no '+str(prof)+' in file '+str(fileEN4)+' is not in the domain, do not keep')
                    i0,j0=-1,-1
                else:
                    print('Profile no '+str(prof)+' in file '+str(fileEN4)+' is in the domain, go proceed')
                    i0=int(last_line.split()[1])
                    j0=int(last_line.split()[2])
    else:
        print('Profile no '+str(prof)+' in file '+str(fileEN4)+' is not in the domain, do not keep')
        i0,j0=-1,-1
       
    return i0,j0


# Check if profile location is on land
def check_prof_in_ocean(i0,j0):
    print('check if profile is in the ocean : ')
    dsN=xr.open_dataset(meshfile)
    tmaskN=dsN.tmask
    if tmaskN[0,0,int(j0)-1,int(i0)-1] == 1:
        check=0
    else:
        check=1
    return check


# Check depth of profile
def check_prof_depth(fileEN4,prof):
    print('check if profile has a good depth : ')
    #open the profile file and read the infos on pressure and latitude
    tfileEN4=diren4+fileEN4
    dsen4=xr.open_dataset(tfileEN4)
    presen4=dsen4['PRES_ADJUSTED'][prof].values
    laten4=dsen4['LATITUDE'][prof].values
    #convert pressure to depth
    depthen4=seawater.dpth(presen4,laten4)
    #look for the last level
    indzprof=np.min(np.where(np.isnan(depthen4)==True))
    dmax=depthen4[indzprof-1]
    print('profile max depth is '+str(dmax)+' m')
    if dmax >= depthmin:
        check=0
    else:
        check=1
    return check


# Check if there are enough model profiles around the obs profile
def map_profile_from_jsonfile(lon,lat,radius,lonmin, lonmax, latmin, latmax,fileEN4,prof):
    # Produce a map with the profile location and the area in which the model profiles are sampled
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.patches as mpatches
    fig=plt.figure(figsize=(20,15))
    ax = plt.subplot(111,projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_extent((lonmin, lonmax, latmin, latmax))
    ax.coastlines(resolution="10m")
    gl = ax.gridlines(draw_labels=True, linestyle=':', color='black',
                      alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.tick_params('both',labelsize=22)
    ax.add_patch(mpatches.Circle(xy=[lon,lat], radius=radius, color='green', alpha=0.3, transform=ccrs.PlateCarree(), zorder=30))
    plt.scatter(lon,lat, c='g', linewidth='0', s=18);
    plt.savefig(plotdir+'/debug_map_'+str(prof)+'_'+str(fileEN4)+'.png')

def check_number_profile(fileEN4,prof,i0,j0,depthmin,dmap=0):
    print('check if there are enough model profiles : ')
    #open mask file and read lat, lon and mask
    ds=xr.open_dataset(meshfile)
    gdpts=np.int(np.round(radius_max*60))
    lat=ds.nav_lat[j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    lon=ds.nav_lon[j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    depth=ds.gdept_1d[0]
    tmask=ds.tmask[0,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
    # Stack the variables
    lon_stacked = lon.stack(profile=('x', 'y'))
    lat_stacked = lat.stack(profile=('x', 'y'))
    mask_stacked = tmask.stack(profile=('x', 'y'))
    #Get the depth at every grid point
    d,ly,lx=tmask.shape
    depthmod2d=np.zeros([lx,ly])
    for j in np.arange(ly):
        for i in np.arange(lx):
            depthmod2d[j,i]=depth[np.min(np.where(tmask[:,j,i].values<1))].values
    xr_depthmod2d=xr.DataArray(depthmod2d, dims=("y", "x"))    
    depth_stacked = xr_depthmod2d.stack(profile=('x', 'y'))
    #open the profile file and read the infos on latitude, longitude
    tfileEN4=diren4+fileEN4
    dsen4=xr.open_dataset(tfileEN4)
    laten4=dsen4['LATITUDE'][prof].values
    lonen4=dsen4['LONGITUDE'][prof].values
    #find the profiles filling criteria
    distance_threshold = radius_max
    square_distance_to_observation = (lon_stacked - lonen4)**2 + (lat_stacked-laten4)**2
    is_close_to_observation = (square_distance_to_observation < distance_threshold) & (depth_stacked > depthmin)
    nb_profiles=np.sum(1*is_close_to_observation)*24*(period*2+1)
    print('There is a total of '+str(nb_profiles)+' model oceanic profiles with enough depth')
    if dmap == 1:
        map_profile_from_jsonfile(lonen4,laten4,radius_max,lonmin, lonmax, latmin, latmax)
    if nb_profiles > number_of_model_profiles:
        check=0
    else:
        check=1
    return check


# Select one profile and save its characteristics in the json file
def make_dict(fileEN4,prof,i0,j0,dictyml):
    print('Adding profile to the dict')
    #open the profile file and read the infos on latitude, longitude
    tfileEN4=diren4+fileEN4
    dsen4=xr.open_dataset(tfileEN4)
    laten4=dsen4['LATITUDE'][prof].values
    lonen4=dsen4['LONGITUDE'][prof].values
    dayen4=dsen4['JULD'][prof].values
    if dictyml:
        up={'Profile_'+str(fileEN4)+'_'+str(prof):{'reference':'Profile_'+str(fileEN4)+'_'+str(prof),
                                                        'file':fileEN4,'profile no':int(prof),
                                                        'latitude':float(laten4),
                                                        'longitude':float(lonen4),
                                                        'date':str(dayen4),
                                                        'i0':int(i0),'j0':int(j0)}}
        dictyml.update(up)
    else:
        dictyml={'Profile_'+str(fileEN4)+'_'+str(prof):{'reference':'Profile_'+str(fileEN4)+'_'+str(prof),
                                                        'file':fileEN4,'profile no':int(prof),
                                                        'latitude':float(laten4),
                                                        'longitude':float(lonen4),
                                                        'date':str(dayen4),
                                                        'i0':int(i0),'j0':int(j0)}}
    return dictyml

# Make the selection of profiles
def selection(ymin,mmin,dmin,ymax,mmax,dmax,period,depthmin,dmap):
    # determining the dates
    datemin=datetime.date(ymin,mmin,dmin)
    datemax=datetime.date(ymax,mmax,dmax)
    # list of days between datemin and datemax
    def date_range(start, end, period):
        r = (end+datetime.timedelta(days=1)-start).days
        return [start+datetime.timedelta(days=i) for i in range(period,r-period)]
    dateList = date_range(datemin, datemax,period) 
    for date in dateList:
        print(date)
    # List of all EN4 files to search in
    list_filesEN4=[]
    for date in dateList:
        year=date.year
        month=date.month
        day=date.day
        mm="{:02d}".format(month) #month on 2 digits
        dd="{:02d}".format(day) # day on 2 digits
        list_filesEN4.append(str(year)+str(mm)+str(dd)+'_prof.nc')
    print(list_filesEN4)
    # loop on all the files
    jdict={}
    for f in range(len(list_filesEN4)):
        fileEN4=list_filesEN4[f]
        ds=xr.open_dataset(diren4+fileEN4)
        nprof=len(ds.N_PROF)
        for prof in range(nprof):
            i0,j0=loc(fileEN4,prof)
            if (i0,j0) == (-1,-1):
                print('profile is not in the domain at all')
                continue
            check=check_prof_in_ocean(i0,j0)
            if check == 1:
                print('no, profile is on the land')
                continue
            print('yes, profile is in the ocean')
            check=check_prof_depth(fileEN4,prof)
            if check == 1:
                print('no, profile is not deep enough')
                continue
            print('yes, profile is deep enough')
            check_number_profile(fileEN4,prof,i0,j0)
            print(check)
            if check == 1:
                print('no, there are not enough model profiles')
                continue
            print('yes, there are enough model profiles')
            jdict=make_dict(fileEN4,prof,i0,j0,jdict)
            # write dict in a json file           
            # name of the json file in which selection of profiles informations will be stored
            jsonfile=jsondir+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'
            with io.open(str(jsonfile), 'a+', encoding='utf8') as outfile:
                outfile.write(str(json.dumps(jdict, sort_keys=True,indent=4, separators=(',', ': '))))
                        
# Plot the locations of all profiles

def plot_profiles_EN4():
    datemin=datetime.date(ymin,mmin,dmin)
    datemax=datetime.date(ymax,mmax,dmax)
    jsonfile=sondir+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'
    sourcefile=open(str(jsonfile),'rU',encoding='utf-8')
    infos=json.load(sourcefile)
    nb_profilesEN4=len(infos)
    list_profiles=infos.keys()

    all_lat=np.zeros((nb_profilesEN4))
    all_lon=np.zeros((nb_profilesEN4))

    for prof in np.arange(nb_profilesEN4):
        reference =  str(list(list_profiles)[prof])
        lat_prof = infos[reference]['latitude']
        lon_prof = infos[reference]['longitude']
        all_lat[prof]=lat_prof
        all_lon[prof]=lon_prof
                        
    # location and name of the maskfile of the model configuration
    ds=xr.open_dataset(meshfile)
    lat=ds.nav_lat
    lon=ds.nav_lon
    tmask=ds.tmask
    dsb=xr.open_dataset(batfile)
    bathy=dsb.Bathymetry
    bathy_mask=np.ma.masked_where(tmask[0,0]==0.,bathy)
    latmin,latmax,lonmin,lonmax=(lat.min(),lat.max(),lon.min(),lon.max())
    datemin=datetime.date(ymin,mmin,dmin)
    datemax=datetime.date(ymax,mmax,dmax)


    fig, axs = plt.subplots(1,2, figsize=(15, 7.5), gridspec_kw={'width_ratios': [4, 1]}, subplot_kw={'projection': ccrs.PlateCarree()})
    axs = axs.ravel()
    axs[0].set_extent((lonmin, lonmax, latmin, latmax))
    pcolor=axs[0].pcolormesh(lon,lat,bathy_mask,transform=ccrs.PlateCarree(),
                             cmap=cmocean.cm.deep,vmin=0,vmax=4000)
    axs[0].coastlines(resolution="10m")
    gl = axs[0].gridlines(draw_labels=True, linestyle=':', color='black',
                          alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    axs[0].tick_params('both',labelsize=22)

    cbar = plt.colorbar(pcolor,orientation='vertical',shrink=0.75,label='m',ax=axs[0])
    axs[0].scatter(all_lon, all_lat, c='r', linewidth='0', s=18);
    axs[0].set_title('There are '+str(len(all_lon))+' EN4 profiles', size=20);

    textstr = '\n'.join((
                ' simulation = MEDWEST60-BLBT02',
                ' dates = '+str(datemin)+' '+str(datemax),
                ' radius max = '+str(radius_max)+'°',
                ' period = '+str(period)+'d',
                ' depth min = '+str(depthmin)+'m',
                ' nb_profiles = '+str(number_of_model_profiles)))        
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axs[1].text(0.05, 0.95, textstr, transform=axs[1].transAxes, fontsize=14,verticalalignment='top', bbox=props)
    axs[1].axis('off')
    fig.tight_layout()

    plt.savefig(plotdir+'/map-profiles-'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.png')

def process_one_profile(prof,infos,dirmod,config,case,meshfile,diren4,radius_max,depthmin,period,number_of_model_profiles,debug_plot=False):
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
        list_filesmod_T=[]
        list_filesmod_S=[]
        for date in dateList:
            year=date.year
            month=date.month
            day=date.day
            mm="{:02d}".format(month) #month on 2 digits
            dd="{:02d}".format(day) # day on 2 digits
            list_filesmod_T.append(dirmod+'/'+str(config)+'-'+str(case)+'_y'+str(year)+'m'+str(mm)+'d'+str(dd)+'.1h_gridT.nc')
            list_filesmod_S.append(dirmod+'/'+str(config)+'-'+str(case)+'_y'+str(year)+'m'+str(mm)+'d'+str(dd)+'.1h_gridS.nc')
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
        tfileEN4=diren4+file_prof
        dsen4=xr.open_dataset(tfileEN4)
        laten4=dsen4['LATITUDE'][prof_prof]
        lonen4=dsen4['LONGITUDE'][prof_prof]
        dayen4=dsen4['JULD'][prof_prof]
        tempen4=dsen4['TEMP_ADJUSTED'][prof_prof]
        salten4=dsen4['PSAL_ADJUSTED'][prof_prof]
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

        # Remove the NaN in the profiles
        indtempnan=np.where(np.isnan(obsred_temp)==True)
        if len(indtempnan[0]) > 0:
            obsred_dep = np.delete(obsred_dep, indtempnan[0])    
            obsred_salt = np.delete(obsred_salt, indtempnan[0])    
            obsred_temp = np.delete(obsred_temp, indtempnan[0])    
        indsaltnan=np.where(np.isnan(obsred_salt)==True)
        if len(indsaltnan[0]) > 0:
            obsred_dep = np.delete(obsred_dep, indsaltnan[0])    
            obsred_salt = np.delete(obsred_salt, indsaltnan[0])    
            obsred_temp = np.delete(obsred_temp, indsaltnan[0])    

        # Find the model profiles
        lon_stacked = lonmod.stack(profile=('x', 'y'))
        lat_stacked = latmod.stack(profile=('x', 'y'))
        mask_stacked = maskmod0.stack(profile=('x', 'y'))
        xr_depthmod2d=xr.DataArray(depthmod2d, dims=('y','x'))
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
            plt.savefig(dirplot+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'_prof'+str(prof)+'.png')

        # Write netcdf file
        match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
        debut_ref=match.group(1)
        fin_ref=match.group(2)
        dirname=diren4+'profiles_files/'+str(config)+'-'+str(case)
        if not os.path.exists(dirname):
            os.makedirs(dirname)    

        outname=dirname+'/'+str(debut_ref)+str(fin_ref)+'_'+str(config)+'-'+str(case)+'_'+str(depthmin)+'m_TS.nc'
        print('output file is '+outname)
        dsout=Dataset(outname,'w')

        today=date.today()
        dsout.description = 'This file contains one profile of temperature and salinity from EN4 dataset and the mean and 10 and 90 percentile of '+str(config)+'-'+str(case)+' data within a '+str(radius_max)+'deg circle around the location of the profile and '+str(period)+' days before and after it has been sampled. This file has been created '+str(today.day)+'/'+str(today.month)+'/'+str(today.year)

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


def process_profiles(ymin,mmin,dmin,ymax,mmax,dmax,config,case,depthmin,radius_max,period,number_of_model_profiles,jsondir,diren4,dirmod,meshfile):
    # Loop over all the profiles listed in the json file                   
    datemin=datetime.date(ymin,mmin,dmin)
    datemax=datetime.date(ymax,mmax,dmax)
    jsonfile=jsondir+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'

    # Read the jsonfile
    sourcefile=open(jsonfile,'rU')
    infos=json.load(sourcefile)
    nb_profilesEN4=len(infos)

    # Loop on the number of profiles
    print("Nb de profiles : "+str(nb_profilesEN4))
    print(time.strftime('%d/%m/%y %H:%M',time.localtime()))

    for prof in range(nb_profilesEN4):
        list_profiles = infos.keys()
        reference =  str(list(list_profiles)[prof])
        match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
        debut_ref=match.group(1)
        fin_ref=match.group(2)
        dirname=diren4+'profiles_files/'+str(config)+'-'+str(case)
        outname=dirname+'/'+str(debut_ref)+str(fin_ref)+'_'+str(config)+'-'+str(case)+'_'+str(depthmin)+'m_TS.nc'
        if not os.path.exists(outname):
            process_one_profile(prof,infos,dirmod,config,case,meshfile,diren4,radius_max,depthmin,period,number_of_model_profiles)
                        
    print(time.strftime('%d/%m/%y %H:%M',time.localtime()))






                
                        
                        


