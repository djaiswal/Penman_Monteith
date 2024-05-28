#  NOTE : This code creates a maskfile to mask the NetCDF files appropriately and retains only the required location.
#         This makes use of a mask file which has value as 1 inside the required locations and value
#         as 0 outside the required locations. This mask file is created using the codes mask0and2flipping.py
#         and mask3652.py

#         Select the absolute filepath to the directory containing the unmasked files from selection window.
#         The second input required is the absolute filepath to the output direcory where the masked files has
#         to be saved.This has to be entered in the terminal.
#         This code also requires the absolute path to the mask file.  
 
import os as os
import netCDF4 as nc
import numpy as np
import sys
import tkinter as tk
from tkinter import filedialog

# Directory containing NetCDF files
directory_path = ""
def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(text=f"Selected Directory with ETo files : {directory_path}")

def browse_maskfile():
        global mask_filepath
        mask_filepath = filedialog.askopenfilename(filetypes=[("NetCDF file", "*.nc"), ("All files", "*.*")])
        label_1.config(text=f"Path to shapefile : {shapefile_path}")

def finish():
    root.destroy()

# Create the main window
root = tk.Tk()
root.title("Directory Browser to mask files in the directory")
root.geometry("800x200")  # Set the window size to 800x200
directory_button = tk.Button(root, text="Browse Directory for masking", command=browse_directory)
directory_button.pack(pady=10)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=20)
finish_button = tk.Button(root, text="Finish selecting directory", command=finish)
finish_button.pack(pady=20)
button_1 = tk.Button(root, text="Browse Mask file", command=browse_maskfile)
button_1.pack(pady=5)
label_1 = tk.Label(root, text="Mask file : ")
label_1.pack(pady=10)
root.mainloop()

output_directory = input("Enter the absolute path to the directory where the masked files are to be stored : ")
# Eg : output_directory = r"E:\ISIMIP Climate Data\et0_masked"  

mask_data = nc.Dataset(mask_filepath, 'r')                                                                                                                                                                                    

for filename in os.listdir(directory_path):
    if filename.endswith('.nc'):
        input_file_path = os.path.join(directory_path, filename)
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
