# All the steps in producing plots of comparison of temp ans salt profiles between ARGO profiles and MEDWEST60 simulations


## Selection of the profiles for the region and date
  - notebook : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-07-21-AA-debug-select-profiles-EN4-region-MED-1month-on-lgge1994.ipynb
  - map of locations of the profiles :
![plot](https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/profiles_MEDWEST60_y2010m02d05-y2010m03d06.png)  
  - problem : get rid of the profiles in Atlantic or near Corsica (too close to the boundary), check if 1° circle in ocean or land in MEDWEST60 mask
  - also : very few profiles ~30 for the region and dates, compare with profiles from other years with same date ? (climatological ?) 
  
=> diag on eNATL60MED 1year of data first !

## Colocalisation of profiles with eNATL60

  - sosie on jean-zay
    - in /gpfswork/rech/egi/rote001/git : ```git clone https://github.com/brodeau/sosie.git```
    - in .bashrc : 
    ```module load intel-all/19.0.4
       module load hdf5/1.8.21-mpi
       module load netcdf/4.7.2-mpi
       module load netcdf-fortran/4.5.2-mpi
       module load netcdf-cxx4/4.3.1-mpi
    ```
    - in /gpfswork/rech/egi/rote001/git/sosie : ``` ln -sf macro/make.macro_ifort_JEAN-ZAY make.macro; make```
  - test ij_from_lon_lat.x with dummy point in MED 
    - script : [test_coloc_sosie.ksh](https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/test_coloc_sosie.ksh)
    - input : [test_coloc.txt]( https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/test_coloc.txt) (no blank line !) and file /gpfsstore/rech/egi/commun/MEDWEST60/
    - output : [ij_found.out](https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/ij_found.out)
