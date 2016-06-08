

input_files = []
# Read netcdf commands
#input_file = '/work/home/sag527/v10/geos5_05x0666_fullchem_ch/HEMCO_diagnostic_copy/HEMCO_diagnostics.200804010200.nc'
input_files.append( 'simple_xy.nc')
#input_files.append( '/work/home/bn506/temp_shani_scripts/planeflight.nc')
#input_files.append( '/work/home/bn506/temp_shani_scripts/shani_planeflight.nc')
##input_file = '/work/home/bn506/Rate_Sensetivity/Batch_runs/CO_OH_p3/results.nc'
#input_files.append('/work/home/bn506/Rate_Sensetivity/Batch_runs/OH_NO2_m_p1/results.nc')
#
#
#input_files.append('/work/home/bn506/Rate_Sensetivity/Monty-Carlo/present/monty_carlo_021/results.nc')

#input_files.append('/work/home/bn506/Rate_Sensetivity/Batch_runs/top_3_p0/ctm.nc')
#input_files.append('/work/home/bn506/planeflight/planeflight.nc')
#input_files.append('/work/home/bn506/test/default_output/ctm.nc')


variable_name = 'latitude'
#variable_name = 'IJ_AVG_S__H2O2'



import netCDF4
from netCDF4 import Dataset
for input_file in input_files:

   print "#####################################################################"
   print input_file
   netCDF_data = Dataset(input_file, mode='r')

   
   print "netCDF data = \n"
   for item in netCDF_data.variables:
       print item
   
   
   
   #variable_data = netCDF_data.variables[variable_name]
   
   print "\n Variable data = \n"
#   print variable_data
   
   print '\n Data = \n'
#   print variable_data[:]
