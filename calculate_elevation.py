# NOTE : This code is used to replicate the Zenith for 3652 timesteps (as we are calculating ETo for 10 years)
#        from the Zenith file for one timestep. This is saved as elevation_data_for_3652_days.nc in the selected directory. 

#        Select the Zenith file with data for one time step.
#        Enter the number of timesteps (3653 if time period is 2091-2100, else 3652 days).
#        Select the directory where the Zenith file for 3652 timesteps has to be saved. 

import numpy as np
import netCDF4 as nc
import tkinter as tk
from tkinter import filedialog
import os

# This part of the code is used to create the interface to select the files and directories
time_size = None  # Initialize

def browse_Zenithfile():
    global Zenith_file_path
    Zenith_file_path = filedialog.askopenfilename(filetypes=[("NetCDF file", "*.nc")])
    label_1.config(text=f"Path to file with data of elevation for one time step : {Zenith_file_path}")

def set_time_size():
    global time_size
    try:
        time_size = int(time_entry.get())
        time_label.config(text=f"Number of timesteps (days): {time_size}")
    except ValueError:
        time_label.config(text="Please enter a valid integer.")

def finish():
    root.destroy()

root = tk.Tk()
root.title("File Browser to select a file with data of elevation")
root.geometry("800x250")
button_1 = tk.Button(root, text="Browse Zenith file", command=browse_Zenithfile)
button_1.pack(pady=5)
label_1 = tk.Label(root, text="Zenith file : ")
label_1.pack(pady=10)

# Entry for time_size
time_label = tk.Label(root, text="Number of timesteps (3653 if time period is 2091-2100, else 3652 days): ")
time_label.pack(pady=5)
time_entry = tk.Entry(root)
time_entry.pack(pady=5)
set_time_button = tk.Button(root, text="Set timesteps", command=set_time_size)
set_time_button.pack(pady=5)

finish_button = tk.Button(root, text="Finish selecting file and timesteps", command=finish)
finish_button.pack(pady=20)
root.mainloop()

if time_size is None:
    time_size = 3652  # fallback default

directory_path = ""
def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(text=f"Selected Directory to save file with data of elevation for 3652 days: {directory_path}")

root = tk.Tk()
root.title("Directory Browser to save file with data of elevation ")
root.geometry("800x200")  # Set the window size to 800x200
directory_button = tk.Button(root, text="Browse Directory to save file with data of elevation for 3652 days", command=browse_directory)
directory_button.pack(pady=10)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=20)
finish_button = tk.Button(root, text="Finish selecting directory", command=finish)
finish_button.pack(pady=20)
root.mainloop()

# This part of the code does the computation. 

Zenith_file = nc.Dataset(Zenith_file_path, 'r')


latitudes =   Zenith_file.variables['lat'][:]
longitudes =   Zenith_file.variables['lon'][:]
values =  Zenith_file.variables['elevation'][:]

# Define dimensions for new file
lat_size = len(latitudes)
lon_size = len(longitudes)


# Creating the new NetCDF file 
output_file = os.path.join(directory_path ,'elevation_data_for_3652_days.nc')                     
with nc.Dataset(output_file, 'w', format='NETCDF4') as nc:
    # Define the dimensions
    nc.createDimension('lat', lat_size)
    nc.createDimension('lon', lon_size)
    nc.createDimension('time', time_size)

    # Create variables
    latitudes_var = nc.createVariable('lat', 'f4', ('lat',))
    longitudes_var = nc.createVariable('lon', 'f4', ('lon',))
    time = nc.createVariable('time', 'i4', ('time',))
    Zenith_data = nc.createVariable('Zenith_data', 'f4', ('time', 'lat', 'lon'))

    # Set latitude and longitude values
    latitudes_var[:] = latitudes
    longitudes_var[:] = longitudes

    # Set time values assuming one value per day starting from 1
    time[:] = np.arange(1, time_size + 1)

    # Create a data array by replicating the values for each time step
    replicated_values = np.tile(values, (time_size, 1, 1))
    replicated_values = replicated_values.astype('float32')  # Convert to float type
    replicated_values = replicated_values.filled(np.nan) 

    Zenith_data[:, :, :] = replicated_values

    Zenith_file.close()

print(f'NetCDF file created successfully.')
