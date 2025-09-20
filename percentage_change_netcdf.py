# NOTE : This code calculates the percentage change of ETo values and saves it as a netCDF file.
#        For this we need NetCDF files with ETo data of calculated without considering
#        the effect of CO2 and ETo data calculated considering the effect of CO2 of each decade

import netCDF4 as nc
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os

directory_path_var_1 = ''
file_paths = [""] * 6
GCM_list = ['ETo with C02 2021-2030','ETo without C02 2021-2030',
            'ETo with C02 2051-2060', 'ETo without C02 2051-2060',
            'ETo with C02 2091-2100', 'ETo without C02 2091-2100']



# The following section of this module is to create an interface to select the files.

def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {GCM_list[file_index]} : {file_path}")
def browse_directory():
    global directory_path_var_1
    path = filedialog.askdirectory()
    if path:
        directory_path_var_1 =  path
        directory_label.config(text=f"Selected Directory with ETo files : {directory_path_var_1}")

def finish():
    root.destroy()

# Create the main window for selecting the NetCDF files.
root = tk.Tk()
root.title("File Browser")
root.geometry("800x450")  # Set the window size to 800x300
labels = []
buttons = []

for i in range(len(file_paths)):

    button = tk.Button(
        root, text=f"Browse File for {GCM_list[i]}", command=lambda i=i: browse_file(i))
    button.pack(pady=5)
    buttons.append(button)

    label = tk.Label(root, text=f"{GCM_list[i]}: ")
    label.pack(pady=5)
    labels.append(label)
directory_button = tk.Button(root, text="Browse Directory to save percentage of relative change in ETo", command=browse_directory)
directory_button.pack(pady=5)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=5)
# Create a button to finish selecting files and close the window
finish_button = tk.Button(root, text="Finish Selecting Files", command=finish)
finish_button.pack(pady=5)
root.mainloop()

def calculate_percentage(withCO2file, withoutCO2file,directory_path_var_1 ):
    # Open the existing NetCDF files containing shortwave and longwave radiation data    
    print("without co2 file, with CO2 file",withoutCO2file, withCO2file)
    withco2 = nc.Dataset(withCO2file, 'r')
    withoutco2 = nc.Dataset(withoutCO2file ,'r')
    if '2021' in withCO2file:
        start_year = "2021-2030"
    elif '2051' in withCO2file:
        start_year = "2051-2060"
    elif '2091' in withCO2file:
        start_year = "2091-2100"

    if ('gfdl' in withCO2file) or 'GFDL' in withCO2file:
        gcm_name = 'GFDL-ESM4'
    elif ('ipsl' in withCO2file) or 'IPSL' in withCO2file:
        gcm_name = 'IPSL-CM6A-LR'
    elif 'mpi' in withCO2file or 'MPI'in withCO2file:
        gcm_name = 'MPI-ESM1-2-HR'
    elif 'mri' in withCO2file or 'MRI'in withCO2file:
        gcm_name = 'MRI-ESM2-0'
    elif 'ukesm' in withCO2file or 'UKESM'in withCO2file:
        gcm_name = 'UKESM1-0-LL'
    # Get dimensions from the input files
    time_dim = len(withoutco2.dimensions['time'])
    latitude_dim = len(withoutco2.dimensions['lat'])
    longitude_dim = len(withoutco2.dimensions['lon'])

    # Create a new NetCDF file
    filename = gcm_name+'_'+start_year+'_percentage_change.nc'
    output_file = os.path.join(directory_path_var_1, filename)                         #Destination of output file
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
    
    if '2021' in withCO2file:
        time_var.units = 'days since 2021-1-1 00:00:00'
    elif '2051' in withCO2file:
        time_var.units = 'days since 2051-1-1 00:00:00'
    elif '2091' in withCO2file:
        time_var.units = 'days since 2091-1-1 00:00:00'
    latitude_var.units = 'degrees_north'
    longitude_var.units = 'degrees_east'

    # Define the variable for net radiation
    percentage_change_var = dataset.createVariable('percentage_change', np.float32, ('time', 'lat', 'lon',))

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
    latitude_var[:] = withoutco2.variables['lat'][:]
    longitude_var[:] = withoutco2.variables['lon'][:]
    percentage_change_var[:] = percentage_change_data

    # Close the NetCDF file
    withoutco2.close()
    withco2.close()
    dataset.close()

    print(f'NetCDF file "{output_file}" created successfully.')


calculate_percentage(file_paths[0], file_paths[1], directory_path_var_1)
calculate_percentage(file_paths[2], file_paths[3], directory_path_var_1)

