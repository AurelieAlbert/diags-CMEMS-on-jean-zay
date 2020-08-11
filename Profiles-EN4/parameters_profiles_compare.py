# directory in which EN4 profiles are stored
diren4='/gpfswork/rech/egi/rote001/EN4/'
# simulation to be compared and where to find it
config='MEDWEST60'
case='BLBT02'
member=
dirmod='/gpfsstore/rech/egi/commun/MEDWEST60/extracted_eNATL60/allv/'
meshfile='/gpfsstore/rech/egi/commun/MEDWEST60/MEDWEST60-I/mesh_mask.nc'
batfile='/gpfsstore/rech/egi/commun/MEDWEST60/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4'
# period considered for the comparison
ymin=2010;mmin=1;dmin=1
ymax=2010;mmax=4;dmax=30
# depth of the desired comparison profile in m
depthmin=1000
# radius of the circle around the profile location in which we take the modeled profiles, in °  
radius_max=0.25
# period of time around the profile sampling date in which we take the modeled profiles, in days
period=5
# minimum amount of model profiles to be considered to make a significant statistical comparison, for instance in a 1° square and 30-days window we have 2.6 millions modeled profiles, in a 0.5°x10 days 216 000
number_of_model_profiles=100000
# where to store all the plots and json files
plotdir='plots'
jsondir='txt'
# wether we wawnt to plot the area of comparison for each profile
dmap=0
