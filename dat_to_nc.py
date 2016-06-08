import pandas as pd
import numpy as np
import netCDF4
from netCDF4 import Dataset
import glob,sys,os,time
import multiprocessing as mp 
 
if len(sys.argv) == 1: sys.argv.append('1')



ncfile = Dataset('data.nc','w')
   
ncfile.initial_conditions_str = ''.join(tuple(open('Init_cons.dat')))
ncfile.date = time.strftime("processing time :%A %d %B %Y at %H:%M")
#ncfile.description = ic_open[0].strip() 

spec = ncfile.createDimension('spec', None)
time = ncfile.createDimension('time', None)
rate = ncfile.createDimension('rate', None)
    
    
def fix_fort (dat):
    dat = dat.strip().split()
    #print dat   
    dummy = []
    for y in dat:
        try: dummy.append(float(y)) 
        except: 
            y = y.replace('-','E-')
            if y[0] == 'E': y = y[1:]
            dummy.append(float(y)) 
    return dummy
    
    
 
for i in sys.argv[1:]:
    
    spc = tuple(open('Spec_%s.dat'%i))
    rte = tuple(open('Spec_%s.dat'%i))
    
    s_head = [i for i in spc[0].strip().replace(' ','').split('!')]
    r_head = [i for i in rte[0].strip().replace(' ','').split('!')]
   
    spc = mp.Pool(16).map(fix_fort, spc[1:])
    rte = mp.Pool(16).map(fix_fort, rte[1:])
    
    
   
    
    
    group = ncfile.createGroup('group_%s'%i)
    specvar = group.createVariable( 'Spec' , "f8"  ,('time','spec',))
    ratevar = group.createVariable( 'Rate' , "f8"  ,('time','rate',))
    
    
    
    
   
    specvar[:] = spc
    ratevar[:] = rte

    specvar.head = ','.join([''.join(x).strip(' ') for x in s_head])
    ratevar.head = ','.join([''.join(x).strip(' ') for x in r_head])

    
ncfile.close()

