# NOTE : This code calculates the seasonal average of ETo for all NetCDF files in a selected directory 
#        and also plots their seasonal plots
#        The first window asks to select the directory where the NetCDF files are present.
#        Next window asks to select the location where the shapefile is located. 

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs                                                
from datetime import datetime, timedelta
import geopandas as gpd
import os
import glob
import pandas as pd
import tkinter as tk
from tkinter import filedialog

# Directory containing NetCDF files
directory_path = ""
def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(text=f"Selected Directory with ETo files : {directory_path}")

def finish():
    root.destroy()

# Create the main window
root = tk.Tk()
root.title("Directory Browser to find seasonal average")
root.geometry("800x200")  # Set the window size to 800x200
directory_button = tk.Button(root, text="Browse Directory with Eto files", command=browse_directory)
directory_button.pack(pady=10)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=20)
finish_button = tk.Button(root, text="Finish selecting directory", command=finish)
finish_button.pack(pady=20)
root.mainloop()
# Define the pattern to match NetCDF files
file_pattern = '*.nc'

# Create a list of file paths
file_paths = glob.glob(os.path.join(directory_path, file_pattern))

# Iterate through each file
for file_path in file_paths:
    # Open the NetCDF file
    dataset = nc.Dataset(file_path)

    time_var = dataset.variables['time']
    eto_var = dataset.variables['evapotranspiration'][:]

    # Extract year and month information from the time variable
    if '2021' in file_path:
        base_date = datetime.strptime('2021-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
    elif '2091' in file_path:
        base_date = datetime.strptime('2091-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
    dates = [base_date + timedelta(days=int(days)) for days in time_var[:]]
    years = np.array([date.year for date in dates])
    months = np.array([date.month for date in dates])

    # Define the time periods for each season
    winter_months = [1, 2]
    pre_monsoon_months = [3, 4, 5]
    monsoon_months = [6, 7, 8, 9]
    post_monsoon_months = [10, 11, 12]

    # Create masks for each season
    winter_mask = np.isin(months, winter_months)
    pre_monsoon_mask = np.isin(months, pre_monsoon_months)
    monsoon_mask = np.isin(months, monsoon_months)
    post_monsoon_mask = np.isin(months, post_monsoon_months)

    # Calculate mean evapotranspiration for each season over the 10 years
    mean_winter_eto = np.nanmean(eto_var[winter_mask, :, :], axis=0)
    mean_pre_monsoon_eto = np.nanmean(eto_var[pre_monsoon_mask, :, :], axis=0)
    mean_monsoon_eto = np.nanmean(eto_var[monsoon_mask, :, :], axis=0)
    mean_post_monsoon_eto = np.nanmean(eto_var[post_monsoon_mask, :, :], axis=0)

    mean_timeandspace_winter_eto = np.nanmean(eto_var[winter_mask, :, :], axis=(0,1,2))
    mean_timeandspace_pre_monsoon_eto = np.nanmean(eto_var[pre_monsoon_mask, :, :], axis=(0,1,2))
    mean_timeandspace_monsoon_eto = np.nanmean(eto_var[monsoon_mask, :, :], axis=(0,1,2))
    mean_timeandspace_post_monsoon_eto = np.nanmean(eto_var[post_monsoon_mask, :, :], axis=(0,1,2))

    print("\n\n")
    print("Seasonal mean of ETo over time and space")
    print("\n")
    if 'gfdl' in file_path:
        print("GCM : GFDL-ESM04")
    elif 'ipsl' in file_path:
        print("GCM : IPSL-CM6A-LR")
    elif 'mpi' in file_path:
        print("GCM : MPI-ESM1-2-HR")
    elif 'mri' in file_path:
        print("GCM : MRI-ESM2-0")
    elif 'ukesm' in file_path:
        print("GCM : UKESM1-0-LL")   

    print("Winter : ", mean_timeandspace_winter_eto)
    print("Pre-monsoon : ", mean_timeandspace_pre_monsoon_eto)
    print("Monsoon : ",mean_timeandspace_monsoon_eto)
    print("Post-monsoon : ", mean_timeandspace_post_monsoon_eto )
        
    mean_et=[]
    mean_et.append(mean_winter_eto)
    mean_et.append(mean_pre_monsoon_eto)
    mean_et.append(mean_monsoon_eto)
    mean_et.append(mean_post_monsoon_eto)

       
    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    time = dataset.variables['time'][:]

    # shapefile_path = r"O:/ISIMIP Climate Data and final results/SHapefiles and netcdf file_02_12_2023/india_state_boundary_projected.shp"
    shapefile_path = ""
    def browse_shapefile():
        global shapefile_path
        shapefile_path = filedialog.askopenfilename(filetypes=[("Shapefiles", "*.shp"), ("All files", "*.*")])
        label_1.config(text=f"Path to shapefile : {shapefile_path}")
    root = tk.Tk()
    root.title("File Browser")
    root.geometry("800x200") 
    button_1 = tk.Button(root, text="Browse shapefile", command=browse_shapefile)
    button_1.pack(pady=5)
    label_1 = tk.Label(root, text="Shapefile : ")
    label_1.pack(pady=10)
    finish_button = tk.Button(root, text="Finish", command=finish)
    finish_button.pack(pady=20)
    root.mainloop()

    states = gpd.read_file(shapefile_path)
    states = states.to_crs(ccrs.PlateCarree().proj4_init)

    max_value =max(np.nanmax(mean_et[0]), np.nanmax(mean_et[1]), np.nanmax(mean_et[2]), np.nanmax(mean_et[3]))
    fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(16, 10))
    plt.subplots_adjust(left=0.275, right=0.725, top=0.9, bottom=0.075, wspace=0.02, hspace=0.2)
    
    for i, ax in enumerate(axes.flat):
        plot = ax.contourf(lon, lat, mean_et[i], cmap='pink_r', transform=ccrs.PlateCarree(), levels=np.linspace(0, max_value, 20))
        ax.set_extent([67.75, 97.75, 7.25, 37.25], crs=ccrs.PlateCarree()) 
        states.boundary.plot(ax=ax, linewidth=0.25, edgecolor='black')
        titles = ['mean_winter_eto','mean_pre_monsoon_eto', 'mean_monsoon_eto', 'mean_post_monsoon_eto']
        ax.set_title(titles[i])
        

    # Create colorbar
    cax0 = fig.add_axes([0.735, 0.25, 0.02, 0.5])  # Position of colorbar for subplot 1
    cbar0 = plt.colorbar(plot, cax=cax0, orientation='vertical', label='Evapotranspiration (in mm/day)', extend='both')

    plt.suptitle(file_path[file_path.find('1_')+2:file_path.find('.nc')].format(max_value),fontsize=18)
    plt.show()
    dataset.close()
