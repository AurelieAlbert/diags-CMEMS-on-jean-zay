# All the steps in producing plots of comparison of temp and salt profiles between ARGO profiles and MEDWEST60 simulations

## Bugs

  - some strange values in the EN4 profiles => take TEMP_ADJUSTED instead of TEMP, get rid of NAN
  - in the script that I want to launch from bash to make it parallelized, the use of xarray when loading multiple files makes the code crash because it cannot launch a new thread (only one proc) + with Dataset it is not possible to do a MFDataset with netcdf4 files => read with a loop on the files
  - in plot final, the mean of profiles was not computed but the last one charged was display (facepalm ...)
  - in processing the profiles, when transforming depth of last level in the model in a data array in order to stack it, x and y were not in the right order so the stack is not the same than for the other varaibles and wrong profiles were taken into account in the mean and percent
  
  
## The process and plot scripts

  - processing the profiles in a json file : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-08-07-AA-process-profiles.ipynb
  - the final plot : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-08-07-AA-final-plot.ipynb


## Redo the selection script

  - a generic script for all the checks and selection of the fitting profiles : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-08-03-AA-selection-script-with-options.ipynb
  - a generic script to produce the map of the selected profiles : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-08-05-AA-map-selected-profiles.ipynb
  - a papermill script to execute those scripts with a set of parameters : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-08-05-AA-execute-map-selected-profiles.ipynb
  - results are a jsonfile in the txt directory and a png in the plots directory
  
## Some preliminary tests to be done before processing the profile

  - check if profile is in the ocean, regarding to the MEDWEST60 mask (if profile in the Atlantic, dismissed)
  - check if profile is taken in the right period : we need X days of model outputs before and after, so profile taken in 01/01/2010 dismissed (put it in the selection of profiles ?)
  - check if profile is too close to the coast or islands : in this case not enough model points but we can think of taking more profiles farther to compensate
  - check if the profile has a good depth : it could be an option, for instance we want a 1000m depth profile so we take all the EN4 profiles that go below 1000m to be sure
      - I tried to install the seawater package via conda-pack because I couldn't append the conda list when load module load conda, kernel is installed ok but the package is not available in notebook (in ipython is fine)
      - install miniconda, create seawater environment from profiles.txt (list of packages)
 
  
-> all these steps in the selection processs, with options :  duration and radius of the model profiles considered, depth

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
    
  - test of ij_from_lon_lat in a notebook : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-07-22-AA-test-call-ij-from-lon-lat-from-sosie-in-notebook.ipynb, test with one profile ok, when all profiles are processed, if one point is not in the domain, it fails
  - localisation on the map with radius : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-07-23-AA-dev-test-process-one-profile.ipynb
  
## Selection of the profiles for the region and date
  - notebook : https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/2020-07-21-AA-debug-select-profiles-EN4-region-MED-1month-on-lgge1994.ipynb
  - map of locations of the 39 profiles :
![plot](https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/profiles_MEDWEST60_y2010m02d05-y2010m03d06.png)  
  - problem : get rid of the profiles in Atlantic or near Corsica (too close to the boundary), check if 1° circle in ocean or land in MEDWEST60 mask
  - also : very few profiles ~30 for the region and dates, compare with profiles from other years with same date ? (climatological ?) 
  
=> diag on eNATL60MED 1year of data first : 307 profiles from 01/01/2010 to 30/09/2010 ![plot](https://github.com/AurelieAlbert/diags-CMEMS-on-jean-zay/blob/master/Profiles-EN4/profiles_MEDWEST60_y2010m01d01-y2010m09d30.png)

   
