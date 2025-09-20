
# NOTE : This code calculates the relative change in ETo under the scenarios
#        when the effect of CO2 is considered and the effect of CO2 is not considered.

#        Select the NetCDF file containing the data of ETo calculated without considering the
#        effect of CO2 and the NetCDF file containing the data of ETo calculated considering the
#        effect of CO2 for the years 2021-2030, 2051-2060, and 2091-2100. Also select the directory 
#        where the NetCDF file with the difference in ETo has to be saved.


import netCDF4 as nc
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os

#  This part of the code creates the interface to select the files. 

# Open the existing NetCDF files
file_paths = [""] * 6
file_list = ['NetCDF file with CO2 2021-2030','NetCDF file without CO2 2021-2030', 
             'NetCDF file with CO2 2051-2060', 'NetCDF file without CO2 2051-2060',
             'NetCDF file  with CO2 2091-2100', 'NetCDF file  without CO2 2091-2100']

def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {file_list[file_index]} : {file_path}")
def finish():
    root.destroy()

directory_path = ""
def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(text=f"Selected Directory to save ETo files : {directory_path}")

root = tk.Tk()
root.title("File Browser")
root.geometry("800x700")  # Set the window size to 800x300

labels = []
buttons = []

for i in range(len(file_paths)):

    button = tk.Button(
        root, text=f"Browse File for {file_list[i]}", command=lambda i=i: browse_file(i))
    button.pack(pady=5)
    buttons.append(button)

    label = tk.Label(root, text=f"{file_list[i]}: ")
    label.pack(pady=5)
    labels.append(label)

# Create a button to finish selecting files and close the window
directory_button = tk.Button(root, text="Select directory to save NetCDF files", command=browse_directory)
directory_button.pack(pady=5)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=5)
finish_button = tk.Button(root, text="Finish Selecting Files", command=finish)
finish_button.pack(pady=5)
root.mainloop()

if ('gfdl' in file_paths[0]) or 'GFDL' in file_paths[0]:
    gcm_name = 'GFDL-ESM4'
elif ('ipsl' in file_paths[0]) or 'IPSL' in file_paths[0]:
    gcm_name = 'IPSL-CM6A-LR'
elif 'mpi' in file_paths[0] or 'MPI'in file_paths[0]:
    gcm_name = 'MPI-ESM1-2-HR'
elif 'mri' in file_paths[0] or 'MRI'in file_paths[0]:
    gcm_name = 'MRI-ESM2-0'
elif 'ukesm' in file_paths[0] or 'UKESM'in file_paths[0]:
    gcm_name = 'UKESM1-0-LL'


if 'pre_monsoon' in file_paths[0]:
    gcm_name = gcm_name + '_pre_monsoon'
elif 'post_monsoon' in file_paths[0]:
    gcm_name = gcm_name + '_post_monsoon'
elif 'winter' in file_paths[0]:
    gcm_name = gcm_name + '_winter'
if 'monsoon' in file_paths[0]:
    gcm_name = gcm_name + '_monsoon'

file_name_2021_30 = gcm_name +'_2021-2030_daily_difference.nc'
file_name_2091_100 = gcm_name +'_2091-2100_daily_difference.nc'
file_name_2051_60 = gcm_name +'_2051-2060_daily_difference.nc'

# Create a new NetCDF file
output_file_2021_30 = os.path.join(directory_path, file_name_2021_30)
output_file_2091_100 = os.path.join(directory_path, file_name_2091_100)
output_file_2051_60 = os.path.join(directory_path, file_name_2051_60)


# This part of the code creates the NetCDF files with the difference in ETo

# Function to create netcdf file
def create_netcdf(output_file, with_CO2_file, without_CO2_file, ):
    dataset = nc.Dataset(output_file, 'w', format='NETCDF4')

    with_CO2 = nc.Dataset(with_CO2_file, 'r')
    without_CO2 = nc.Dataset(without_CO2_file, 'r')
    
    time_dim = len(with_CO2.dimensions['time'])                                          # size of time dimension
    latitude_dim = len(with_CO2.dimensions['lat'])                                        # latitude dimension
    longitude_dim = len(with_CO2.dimensions['lon'])                                      # longitude dimension

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
    elif '2051' in with_CO2_file:
        time_var.units = 'days since 2051-1-1 00:00:00'
    latitude_var.units = 'degrees_north'
    longitude_var.units = 'degrees_east'

    # Define the variable for net radiation
    output_var = dataset.createVariable(
        'difference_in_eto', np.float32, ('time', 'lat', 'lon', ))

    # Define units and attributes for net radiation
    output_var.units = 'mm/day'
    output_var.long_name = 'Difference in ETo with and without the effect of CO2'

    withoutCO2_data = without_CO2.variables['evapotranspiration'][:]
    withCO2_data = with_CO2.variables['evapotranspiration'][:]

    # Calculate the difference in ETo
    difference_data = withoutCO2_data - withCO2_data

    # Assign data to variables
    time_var[:] = without_CO2.variables['time'][:]
    latitude_var[:] = without_CO2.variables['lat'][:]
    longitude_var[:] = without_CO2.variables['lon'][:]
    # output_var[:] = np.nanmean(difference_data, axis =0)
    output_var[:] = difference_data

    # Close the NetCDF file
    without_CO2.close()
    with_CO2.close()
    dataset.close()

    print(f'NetCDF file "{output_file}" created successfully.')
    return

create_netcdf(output_file = output_file_2051_60, with_CO2_file = file_paths[0], without_CO2_file = file_paths[1])
create_netcdf(output_file = output_file_2051_60, with_CO2_file = file_paths[2], without_CO2_file = file_paths[3])
create_netcdf(output_file = output_file_2091_100, with_CO2_file = file_paths[4], without_CO2_file = file_paths[5])
