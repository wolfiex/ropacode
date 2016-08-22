import numpy as np
import pandas as pd
import sys,os,re,multiprocessing,netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#netcdf file path
ncfile = sys.argv[1]


netCDF_data = Dataset(ncfile, mode='r')

for group in netCDF_data.groups:
    print '\n\nReading', group
    #data
    spc = pd.DataFrame(netCDF_data.groups[group].variables['Spec'][:])
    spc.columns = str(netCDF_data.groups[group].variables['Spec'].head).split(',')

    #rte = pd.DataFrame(netCDF_data.groups[group].variables['Rate'][:])
    #rte.columns = str(netCDF_data.groups[group].variables['Rate'].head).split(',')


    spc = spc/spc.M.mean()
    spc.sort_index(axis=1,inplace=True)# arrange alphabetically
    
      
    
    pp = PdfPages('%s.pdf'%group)
    
    
    for i in xrange(0, len(spc.columns), 6):
        spc[spc.columns[i:i+5]].plot(subplots=True)
        plt.tight_layout() 
        plt.ylabel('mix ratio')
        
        #plt.locator_params(axis='y',nbins=2)
        print '%.03f'%(float(i) / float(len(spc.columns)) ) , '% done'  
        plt.savefig(pp, format='pdf')
        plt.close()
        
    pp.close()
    print 'PDF out'
