
# NOTE : This code code calculates the relative change in ETo under the scenarios
#        when the effect of CO2 is considered and the effect of CO2 is not considered.

#        The first input required is the absolute filepath to the NetCDF file containing
#        the data of ETo calculated without considering the effect of CO2.
#        The second input required is the absolute filepath to the NetCDF file containing
#        the data of ETo calculated considering the effect of CO2.
#        The third input required is the absolute filepath to the NetCDF where the the data
#        corresponding to the difference in ETo (when CO2 is not considered and CO2 is considered)
#        has to be stored. (This file need not exist before running this code. It will be created if
#        it does not exist already.)


import netCDF4 as nc
import numpy as np

# Open the existing NetCDF files

without_CO2_file = input(
    "Enter the absolute filepath to the NetCDF file with ETo data calculated without CO2 : ")
# Eg: without_CO2_file = r"E:\ISIMIP Climate Data\eto_masked_changed\masked1_ipsl_et0_2021-2030_without_CO2.nc"
with_CO2_file = input(
    "Enter the absolute filepath to the NetCDF file with ETo data calculated considering CO2 : ")
# Eg: with_CO2_file = r"E:\ISIMIP Climate Data\eto_masked_changed\masked1_ipsl_et0_2021-2030_with_CO2.nc"

without_CO2 = nc.Dataset(without_CO2_file, 'r')
with_CO2 = nc.Dataset(with_CO2_file, 'r')

time_dim = 3652                                          # size of time dimension
latitude_dim = 61                                        # latitude dimension
longitude_dim = 61                                       # longitude dimension

# Create a new NetCDF file
output_file = input(
    "Enter the absolute filepath to the NetCDF file where the data has to be stored : ")
# Eg: output_file = r"E:\ISIMIP Climate Data\eto_masked_changed\masked_ipsl_et0_2021-2030_difference.nc"
dataset = nc.Dataset(output_file, 'w', format='NETCDF4')

# Define dimensions in the NetCDF file
time = dataset.createDimension('time', time_dim)
latitude = dataset.createDimension('lat', latitude_dim)
longitude = dataset.createDimension('lon', longitude_dim)

# Create variables for time, latitude, and longitude
time_var = dataset.createVariable('time', np.float32, ('time',))
latitude_var = dataset.createVariable('lat', np.float32, ('lat',))
longitude_var = dataset.createVariable('lon', np.float32, ('lon',))

# Define units for variables
if '2021' in with_CO2_file:
    time_var.units = 'days since 2021-1-1 00:00:00'
elif '2091' in with_CO2_file:
    time_var.units = 'days since 2091-1-1 00:00:00'
latitude_var.units = 'degrees_north'
longitude_var.units = 'degrees_east'

# Define the variable for net radiation
output_var = dataset.createVariable(
    'difference_in_eto', np.float32, ('time', 'lat', 'lon',))

# Define units and attributes for net radiation
output_var.units = 'mm/day'
output_var.long_name = 'Difference in ETo with and without the effect of CO2'

withoutCO2_data = without_CO2.variables['evapotranspiration'][:]
withCO2_data = with_CO2.variables['evapotranspiration'][:]


# Applying mask by multiplying with the mask file
difference_data = withoutCO2_data - withCO2_data

# Assign data to variables
time_var[:] = np.arange(time_dim)
latitude_var[:] = np.linspace(37.25, 7.25, latitude_dim)
longitude_var[:] = np.linspace(67.75, 97.75, longitude_dim)
output_var[:] = difference_data

# Close the NetCDF file
without_CO2.close()
with_CO2.close()
dataset.close()

print(f'NetCDF file "{output_file}" created successfully.')
