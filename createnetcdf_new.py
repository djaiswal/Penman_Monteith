import netCDF4 as nc
import numpy as np

# Open the existing NetCDF files containing shortwave and longwave radiation data
shortwave_file = nc.Dataset(r"D:\GFDL_ESM4_new\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rsds_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc" ,'r')
longwave_file = nc.Dataset(r"D:\GFDL_ESM4_new\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rlds_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc", 'r')

# Get dimensions from the input files
time_dim = len(shortwave_file.dimensions['time'])
latitude_dim = len(shortwave_file.dimensions['lat'])
longitude_dim = len(shortwave_file.dimensions['lon'])

# Create a new NetCDF file
output_file = 'D:\\folder\net_radiationnew'
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
time_var.units = 'days since 2021-1-1 00:00:00'
latitude_var.units = 'degrees_north'
longitude_var.units = 'degrees_east'

# Define the variable for net radiation
net_radiation_var = dataset.createVariable('net_radiation', np.float32, ('time', 'latitude', 'longitude',))

# Define units and attributes for net radiation
net_radiation_var.units = 'W/m^2'
net_radiation_var.long_name = 'Net Radiation (Shortwave + Longwave)'

# Read shortwave and longwave radiation data from the existing files
shortwave_data = shortwave_file.variables['rsds'][:]
longwave_data = longwave_file.variables['rlds'][:] 

# Calculate net radiation by summing shortwave and longwave
net_radiation_data = shortwave_data + longwave_data

# Assign data to variables
time_var[:] = np.arange(time_dim)
latitude_var[:] = shortwave_file.variables['lat'][:]
longitude_var[:] = shortwave_file.variables['lon'][:]
net_radiation_var[:] = net_radiation_data

# Close the NetCDF file
shortwave_file.close()
longwave_file.close()
dataset.close()

print(f'NetCDF file "{output_file}" created successfully.')
