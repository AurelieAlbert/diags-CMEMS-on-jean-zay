#!/usr/bin/env python

import comp_en4_profiles_noxarray as cep
import parameters_profiles_compare_100m as param
import datetime
import json
import re
import os


diren4=param.diren4
config=param.config
case=param.case
member=param.member
dirmod=param.dirmod
meshfile=param.meshfile
batfile=param.batfile
ymin=param.ymin;mmin=param.mmin;dmin=param.dmin
ymax=param.ymax;mmax=param.mmax;dmax=param.dmax
depthmin=param.depthmin
radius_max=param.radius_max
period=param.period
number_of_model_profiles=param.number_of_model_profiles
plotdir=param.plotdir
jsondir=param.jsondir
dmap=param.dmap

datemin=datetime.date(ymin,mmin,dmin)
datemax=datetime.date(ymax,mmax,dmax)
jsonfile=jsondir+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'

sourcefile=open(jsonfile,'rU')
infos=json.load(sourcefile)

prof=0

list_profiles = infos.keys()
print(len(list_profiles))

reference =  str(list(list_profiles)[prof])
match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
debut_ref=match.group(1)
fin_ref=match.group(2)
dirname=diren4+'profiles_files/'+str(config)+'-'+str(case)
outname=dirname+str(debut_ref)+str(fin_ref)+'_'+str(config)+'-'+str(case)+'_'+str(depthmin)+'m_TS.nc'
print(outname)
if not os.path.exists(outname):
    cep.process_one_profile(prof,infos,dirmod,config,case,meshfile,diren4,radius_max,depthmin,period,number_of_model_profiles)
