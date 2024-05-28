# NOTE : This code is used to flip the values to create a proper mask file.
#        It requires the absolute file path to the netCDF file which has value as 0 for required locations
#        and value as 1 for other locations.This file can be obtained by processing the shapefile of 
#        of India in the ArcGIS software. Here the required location is the area in the Indian Mainland.

#        This code requires the absolute file path to the NetCDF file with value as 0 for required locations
#        as the first input.
#        The absolute file path to the the NetCDF file where the new mask file has to be saved is the 
#        second input required.

       
import netCDF4 as nc
import numpy as np
import tkinter as tk
from tkinter import filedialog

def browse_maskfile():
        global mask_file_path
        mask_file_path = filedialog.askopenfilename(filetypes=[("NetCDF file", "*.nc"), ("All files", "*.*")])
        label_1.config(text=f"Path to shapefile : {mask_file_path}")

def finish():
    root.destroy()

root = tk.Tk()
root.title("File Browser to select a maskfile")
root.geometry("800x200")
finish_button = tk.Button(root, text="Finish selecting file", command=finish)
finish_button.pack(pady=20)
button_1 = tk.Button(root, text="Browse Mask file", command=browse_maskfile)
button_1.pack(pady=5)
label_1 = tk.Label(root, text="Mask file : ")
label_1.pack(pady=10)
root.mainloop()

input_mask_file = nc.Dataset(mask_file_path, 'r')

latitude_dim = len(input_mask_file.dimensions['lat'])
longitude_dim = len(input_mask_file.dimensions['lon'])

output_file = input("Enter the absolute filepath to the NetCDF file where the mask file has to be saved : ")
# Eg: output_file = r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1_proper_india_02_12_2023.nc"                        #Destination of output file
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