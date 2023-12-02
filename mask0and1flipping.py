import netCDF4 as nc
import numpy as np

input_mask_file = nc.Dataset(r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1india_02_12_2023.nc", 'r')

latitude_dim = len(input_mask_file.dimensions['lat'])
longitude_dim = len(input_mask_file.dimensions['lon'])

output_file = r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1_proper_india_02_12_2023.nc"                        #Destination of output file
dataset = nc.Dataset(output_file, 'w', format='NETCDF4')

latitude = dataset.createDimension('lat', latitude_dim)
longitude = dataset.createDimension('lon', longitude_dim)

latitude_var = dataset.createVariable('lat', np.float32, ('lat',))
longitude_var = dataset.createVariable('lon', np.float32, ('lon',))
new_var = dataset.createVariable('mask_var', np.float32, ( 'lat', 'lon',))

latitude_var.units = 'degrees_north'
longitude_var.units = 'degrees_east'
new_var.units='mask'

latitude_var[:] = (input_mask_file.variables['lat'][:])
longitude_var[:] = (input_mask_file.variables['lon'][:])
original_array = np.where(input_mask_file.variables['0and1'][:] == 0,1,0)
modified_array = np.where(original_array == 0, np.nan, original_array)

new_var[:] = modified_array
input_mask_file.close()

print(f'NetCDF file "{output_file}" created successfully.')