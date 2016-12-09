'''
A tool to calculate the fluxes from DSMACC
D.Ellis 2016
'''


#functions
global specs,reactants



xlen = lambda x: xrange(len(x))

def getcoef (x): #  '''gets specie coefficients from data'''
    try: return int(re.sub(r'([\.\d]*)\s*\D[\d\D]*', r'\1', x))
    except: return 1  #assume coeff are Z+

def getspec (x): #makes concentration array
    return [[specs.loc[x,re.sub(r'([\.\d\s]*)(\D[\d\D]*)', r'\2', i)] for i in j] for j in reactants]


# make a single large array





import numpy as np
import pandas as pd
import sys,os,re,multiprocessing,netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#netcdf file path
ncfile = sys.argv[1]
ncores = 16


########### read dsmacc data
myfile = 'nheptane_1612091358.nc'


nc = Dataset(myfile,'r')
print nc.date, '\n', nc.description,'\n'
print 'Select Simulation: \n\n'
for i,g in enumerate(nc.groups): print i , ' - ', g
group = tuple(nc.groups)[int(input('Enter Number \n'))]
0
print group, 'took', nc.groups[group].WALL_time, 'seconds to compute.'
specs = pd.DataFrame(nc.groups[group].variables['Spec'][:])
specs.columns = nc.groups[group].variables['Spec'].head.split(',')
rates = pd.DataFrame(nc.groups[group].variables['Rate'][:])
rates.columns = nc.groups[group].variables['Rate'].head.split(',')
print 'Spec and Rate files loaded'



nc.close()
########################
specs['TIME'] = pd.to_datetime(specs.TIME, unit='s')
rates['TIME'] = specs['TIME']





''' selection type here '''
timestep = 1
# designed to run for whole array, although if a selection is to be made, this should be done here

''' 1 remove dummy and non-reactions'''
rates = rates[[r for r in rates.columns[6:] if ('DUMMY' not in r) & ('EMISS' not in r)]]

''' 2 remove species if not present (shrink data) '''
#specs = specs[specs.columns[specs.sum()>=0]]
rates = rates[rates.columns[rates.sum()>=0]]

''' 3 get conversion factor from molcm-3 to mixing ratio'''
M = specs['M'].mean()

''' 4 generate reactants and products list '''
#no nead to clear whitespace as begin.py should take care of that.
rate_head = '\n'+'\n'.join(rates.columns)+'\n'
products = [i.split('+') for i in re.findall(r'-->([A-z0-9+]*)',rate_head)]
reactants = [j.split('+') for j in re.findall(r'\n([A-z0-9+]{1,60})[-->]{0,1}',rate_head)]
if len(reactants) != len(products) : print 'reactants and poducts differing lengths'

''' 5 get reaction coefficients'''
coef_re,coef_pro=[],[]
for i in xlen(reactants):
    coef_re.append(np.vectorize(getcoef)(reactants[i]))

    #coef_pro.append(np.vectorize(getcoef)(products[i])) # unused as flux calucated from reactants
    #products[i] = np.vectorize(getspec)(products[i])

''' 6 trip to only timesteps required '''




''' 7 concentration array '''
conc = np.vectorize(getspec)(xlen(specs))
