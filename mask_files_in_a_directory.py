#  NOTE : This code creates a maskfile to mask the NetCDF files appropriately and retains only the required location.
#         This makes use of a mask file which has value as 1 inside the required locations and value
#         as 0 outside the required locations. This mask file is created using the codes mask0and2flipping.py
#         and mask3652.py

#         This code requires the absolute filepath to the directory containing the unmasked files as 
#         first input. 
#         The second input required is the absolute filepath to the output direcory where the masked files has
#         to be saved.
#         This code also requires the absolute path to the mask file.  
 
import os as os
import netCDF4 as nc
import numpy as np
import sys

input_directory = input("Enter the absolute file path to the directory containing the unmasked NetCDF files : ")
# Eg : input_directory = 'E:\ISIMIP Climate Data\et0_unmasked_7.25-37.25_67.75-97.75'
output_directory = input("Enter the absolute path to the directory where the masked files are to be stored : ")
# Eg : output_directory = r"E:\ISIMIP Climate Data\et0_masked"  
mask_data = input("Enter the absolute filepath to the mask file : ")             # file which has 1 at places inside 
                                                                                 # India and 0 outside India.
# Eg : mask_data= "E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1_proper_lat_lon_adj_india_02_12_2023.nc"
mask_data = nc.Dataset(mask_data, 'r')                                                                                                                                                                                    

for filename in os.listdir(input_directory):
    if filename.endswith('.nc'):
        input_file_path = os.path.join(input_directory, filename)
        output_file_path = os.path.join(output_directory, f'masked_{filename}')

        nc_dataset = nc.Dataset(input_file_path, 'r')


        min_lat, max_lat = 7.25,37.25                                            # minimum and maximum latitudes of India
        min_lon, max_lon = 67.75,97.75                                           # minimum and maximum latituses of India


        lat = nc_dataset.variables['lat'][:]
        lon = nc_dataset.variables['lon'][:]


        lat_indices = np.where((lat >= min_lat) & (lat <= max_lat))[0]
        lon_indices = np.where((lon >= min_lon) & (lon <= max_lon))[0]


        subset_data = nc_dataset.variables[list(nc_dataset.variables.keys())[3]][:, lat_indices, lon_indices] * mask_data.variables['mask_data'][:, :,:]

        # Define dimensions and attributes for output file.
        output_dataset = nc.Dataset(output_file_path, 'w', format='NETCDF4')


        output_dataset.createDimension('lat', len(lat_indices))
        output_dataset.createDimension('lon', len(lon_indices))
        output_dataset.createDimension('time', subset_data.shape[0])


        output_lat = output_dataset.createVariable('lat', 'f4', ('lat',))
        output_lon = output_dataset.createVariable('lon', 'f4', ('lon',))
        output_time = output_dataset.createVariable('time', 'f4', ('time',))
        output_data = output_dataset.createVariable(list(nc_dataset.variables.keys())[3], 'f4', ('time', 'lat', 'lon'))

        if '2021' in filename:
            output_time.units = 'days since 2021-1-1 00:00:00'
        elif '2091' in filename:
            output_time.units = 'days since 2021-1-1 00:00:00'

        output_lat.units = 'degrees_north'
        output_lon.units = 'degrees_east'

        output_lat[:] = lat[lat_indices]
        output_lon[:] = lon[lon_indices]
        output_time[:] = nc_dataset.variables['time'][:]  


        output_data[:, :, :] = subset_data


        nc_dataset.close()
        output_dataset.close()

        print("Subset data has been saved to", output_file_path)
