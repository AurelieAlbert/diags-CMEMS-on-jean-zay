#!/usr/bin/env python
# coding: utf-8

# ## Selection of EN4 profiles to compare to a simulation

# In this script the EN4 database is browsed in order to find the profiles that will be compare to the outputs of one simulation. 
# 
# The comparison will not be a simple colocation of the profile inside the model grid but we want to make a statistical comparison of the observed profile with  a significant number of profiles close to it in the model (close in terms of space and time, for instance in a 0,5° radius around the profile location and 10 days before and after it has been sampled)
# 
# Therefor the selected profiles must have to be relevant according to some criteria :
#   - they must be inside the domain of simulation (the mask file of the configuration file must be provided)  
#   - they must be sampled inside the period of simulation (the period shortened by a certain amount of days must be provided)
#   - they must go as deep as a given depth (according to the desired depth for the comparison profiles)
#   - they must be in a location where there are enough model profiles (for instance if profile is too close to an island or the coast)

# ### Criteria of profiles selection

# In[ ]:


# period of the simulation (year,month,day)
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


# In[ ]:


# location and name of the maskfile of the model configuration
meshfile='/gpfsstore/rech/egi/commun/MEDWEST60/MEDWEST60-I/mesh_mask.nc'


# ### imports of librairies

# In[ ]:


import numpy as np
import dask
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import datetime
import glob as glob
import time
from datetime import date
import io
import json
import seawater
import sys
import os


# ### determining the dates

# In[ ]:


datemin=datetime.date(ymin,mmin,dmin)
datemax=datetime.date(ymax,mmax,dmax)
# list of days between datemin and datemax
def date_range(start, end, period):
    r = (end+datetime.timedelta(days=1)-start).days
    return [start+datetime.timedelta(days=i) for i in range(period,r-period)]
dateList = date_range(datemin, datemax,period) 
for date in dateList:
    print(date)


# ### List of all EN4 files to search in

# In[ ]:


diren4="/gpfswork/rech/egi/rote001/EN4/"
list_filesEN4=[]
for date in dateList:
    year=date.year
    month=date.month
    day=date.day
    mm="{:02d}".format(month) #month on 2 digits
    dd="{:02d}".format(day) # day on 2 digits
    list_filesEN4.append(str(year)+str(mm)+str(dd)+'_prof.nc')
print(list_filesEN4)


# ### Define the criteria for one profile

# #### Localisation of the profile location inside model domain

# In[ ]:


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
            with open('txt/prof0.txt','w') as txt_file:
                txt_file.write('Profile_'+str(fileEN4)+'_'+str(prof)+' '+str(lonen4)+' '+str(laten4))
            if os.path.exists('ij_found.out'):
                get_ipython().system('rm ij_found.out')
            
            get_ipython().system('/gpfswork/rech/egi/rote001/git/sosie/bin/ij_from_lon_lat.x -i $meshfile -p txt/prof0.txt > txt/output')

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
#        print('Profile coordinates are '+str(lonen4)+'°E x '+str(laten4)+'°N')
#        print('Model domain is '+str(latmin.values)+'-'+str(latmax.values)+'°E x '+str(lonmin.values)+'-'+str(lonmax.values)+'°N')
        i0,j0=-1,-1
       
    return i0,j0


# #### Check if profile location is on land

# In[ ]:


def check_prof_in_ocean(i0,j0):
    print('check if profile is in the ocean : ')
    dsN=xr.open_dataset(meshfile)
    tmaskN=dsN.tmask
    if tmaskN[0,0,int(j0)-1,int(i0)-1] == 1:
        check=0
    else:
        check=1
    return check


# #### check depth of profile

# In[ ]:


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


# #### check if there are enough model profiles around the obs profile

# In[ ]:


import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patches as mpatches

def map_profile_from_jsonfile(lon,lat,radius,lonmin, lonmax, latmin, latmax):
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

def check_number_profile(fileEN4,prof,i0,j0,dmap=0):
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
    xr_depthmod2d=xr.DataArray(depthmod2d, dims=("x", "y"))    
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


# #### function to select one profile and save its characteristics in the json file

# In[ ]:


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
    
    


# In[ ]:


get_ipython().run_cell_magic('time', '', "#loop on all the files\n\njdict={}\nfor f in range(len(list_filesEN4)):\n    fileEN4=list_filesEN4[f]\n    ds=xr.open_dataset(diren4+fileEN4)\n    nprof=len(ds.N_PROF)\n    for prof in range(nprof):\n        i0,j0=loc(fileEN4,prof)\n        if (i0,j0) == (-1,-1):\n            print('profile is not in the domain at all')\n            continue\n        check=check_prof_in_ocean(i0,j0)\n        if check == 1:\n            print('no, profile is on the land')\n            continue\n        print('yes, profile is in the ocean')\n        check=check_prof_depth(fileEN4,prof)\n        if check == 1:\n            print('no, profile is not deep enough')\n            continue\n        print('yes, profile is deep enough')\n        check_number_profile(fileEN4,prof,i0,j0)\n        print(check)\n        if check == 1:\n            print('no, there are not enough model profiles')\n            continue\n        print('yes, there are enough model profiles')\n        jdict=make_dict(fileEN4,prof,i0,j0,jdict)\n    \n#write dict in a json file           \n# name of the json file in which selection of profiles informations will be stored\njsonfile='txt/MEDWEST60-BLBT02_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'\nwith io.open(str(jsonfile), 'a+', encoding='utf8') as outfile:\n    outfile.write(str(json.dumps(jdict, sort_keys=True,indent=4, separators=(',', ': '))))")

