import numpy as np
import netCDF4 as nc

mask_file = nc.Dataset(r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1_proper_india_02_12_2023.nc", 'r')


latitudes =   mask_file.variables['lat'][:]
longitudes =   mask_file.variables['lon'][:]
values =  mask_file.variables['mask_var'][:]

# Define dimensions for new file
lat_size = len(latitudes)
lon_size = len(longitudes)
time_size = 3652

# Creating the new NetCDF file 
with nc.Dataset(r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1_proper_india_3652.nc", 'w', format='NETCDF4') as nc:
    # Define the dimensions
    nc.createDimension('lat', lat_size)
    nc.createDimension('lon', lon_size)
    nc.createDimension('time', time_size)

    # Create variables
    latitudes_var = nc.createVariable('lat', 'f4', ('lat',))
    longitudes_var = nc.createVariable('lon', 'f4', ('lon',))
    time = nc.createVariable('time', 'i4', ('time',))
    mask_data = nc.createVariable('mask_data', 'f4', ('time', 'lat', 'lon'))

    # Set latitude and longitude values
    latitudes_var[:] = latitudes
    longitudes_var[:] = longitudes

    # Set time values assuming one value per day starting from 1
    time[:] = np.arange(1, time_size + 1)

    # Create a data array by replicating the values for each time step
    replicated_values = np.tile(values, (time_size, 1, 1))
    mask_data[:, :, :] = replicated_values

    mask_file.close()

print(f'NetCDF file created successfully.')
