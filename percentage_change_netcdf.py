# NOTE : This code calculates the percentage change of ETo values and saves it as a netCDF file.
#        For this we need NetCDF files with ETo data of calculated without considering
#        the effect of CO2 and ETo data of 2021-2030 calculated considering the effect of CO2

import netCDF4 as nc
import numpy as np

# Open the existing NetCDF files containing shortwave and longwave radiation data
withoutco2 = nc.Dataset(r"D:\GFDL_ESM4_new\evapotranspiration_without_co2_2091-2100.nc" ,'r')
withco2 = nc.Dataset(r"D:\GFDL_ESM4_new\evapotranspiration_with_co2_2091-2100.nc", 'r')

# Get dimensions from the input files
time_dim = len(withoutco2.dimensions['time'])
latitude_dim = len(withoutco2.dimensions['latitude'])
longitude_dim = len(withoutco2.dimensions['longitude'])

# Create a new NetCDF file
output_file = 'D:\GFDL_ESM4_new\percentage_change_2091-2100'                         #Destination of output file
dataset = nc.Dataset(output_file, 'w', format='NETCDF4')

# Define dimensions in the NetCDF file
time = dataset.createDimension('time', time_dim)
latitude = dataset.createDimension('latitude', latitude_dim)
longitude = dataset.createDimension('longitude', longitude_dim)

# Create variables for time, latitude, and longitude
time_var = dataset.createVariable('time', np.float32, ('time',))
latitude_var = dataset.createVariable('latitude', np.float32, ('latitude',))
longitude_var = dataset.createVariable('longitude', np.float32, ('longitude',))

# Define units for variables
time_var.units = 'days since 2091-1-1 00:00:00'
latitude_var.units = 'degrees_north'
longitude_var.units = 'degrees_east'

# Define the variable for net radiation
percentage_change_var = dataset.createVariable('percentge_change', np.float32, ('time', 'latitude', 'longitude',))

# Define units and attributes for net radiation
percentage_change_var.units = '%'
percentage_change_var.long_name = 'Difference in percentage'

# Read shortwave and longwave radiation data from the existing files
withoutco2_data = withoutco2.variables['evapotranspiration'][:]
withco2_data = withco2.variables['evapotranspiration'][:] 

# Calculate net radiation by summing shortwave and longwave
percentage_change_data = (withoutco2_data-withco2_data)/(withoutco2_data)*100

# Assign data to variables
time_var[:] = np.arange(time_dim)
latitude_var[:] = withoutco2.variables['latitude'][:]
longitude_var[:] = withoutco2.variables['longitude'][:]
percentage_change_var[:] = percentage_change_data

# Close the NetCDF file
withoutco2.close()
withco2.close()
dataset.close()

print(f'NetCDF file "{output_file}" created successfully.')
