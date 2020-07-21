# All the steps in producing plots of comparison of temp ans salt profiles between ARGO profiles and MEDWEST60 simulations


## Selection of the profiles for the region and date
  - notebook : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-07-21-AA-debug-select-profiles-EN4-region-MED-1month-on-lgge1994.ipynb
  - map of locations of the profiles :
![plot](https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/profiles_MEDWEST60_y2010m02d05-y2010m03d06.png)  
  - problem : get rid of the profiles in Atlantic or near Corsica (too close to the boundary), check if 1° circle in ocean or land in MEDWEST60 mask
  - also : very few profiles ~30 for the region and dates, compare with profiles from other years with same date ? (climatological ?) 
  
=> diag on eNATL60MED 1year of data first !
