# NOTE : This module plots the mean of relative change in ETo for five GCMs for a time period.
#        It is required to select the shapefile and the five different NetCDF files which has the data
#        of difference in daily ETo (difference of ETo calculated with and without considering CO2)
#        for the five GCMs.

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
from netCDF4 import Dataset
import numpy as np
import tkinter as tk
from tkinter import filedialog


# The following section of this module is to create an interface to select the files.
file_paths = [""] * 5
GCM_list = ['GFDL-ESM4', 'IPSL-CM6A-LR',
            'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL']

def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {GCM_list[file_index]} : {file_path}")

def browse_shapefile():
        global shapefile_path
        shapefile_path = filedialog.askopenfilename(filetypes=[("Shapefiles", "*.shp"), ("All files", "*.*")])
        label_1.config(text=f"Path to shapefile : {shapefile_path}")


def finish():
    root.destroy()

root = tk.Tk()
root.title("File Browser")
root.geometry("800x500")  # Set the window size to 800x300

labels = []
buttons = []
button_1 = tk.Button(root, text="Browse shapefile", command=browse_shapefile)
button_1.pack(pady=5)
label_1 = tk.Label(root, text="Shapefile : ")
label_1.pack(pady=10)

for i in range(len(file_paths)):

    button = tk.Button(
        root, text=f"Browse File for {GCM_list[i]}", command=lambda i=i: browse_file(i))
    button.pack(pady=5)
    buttons.append(button)

    label = tk.Label(root, text=f"{GCM_list[i]}: ")
    label.pack(pady=5)
    labels.append(label)

# Create a button to finish selecting files and close the window
finish_button = tk.Button(root, text="Finish Selecting Files", command=finish)
finish_button.pack(pady=5)
root.mainloop()

states = gpd.read_file(shapefile_path)
states = states.to_crs(ccrs.PlateCarree().proj4_init)

gfdl_diff = Dataset(file_paths[0],'r')
ipsl_diff = Dataset(file_paths[1],'r')
mpi_diff = Dataset(file_paths[2],'r')
mri_diff = Dataset(file_paths[3],'r')
ukesm_diff = Dataset(file_paths[4],'r')

gfdl_diff_data = gfdl_diff.variables['difference_in_eto'][:]
ipsl_diff_data = ipsl_diff.variables['difference_in_eto'][:]
mpi_diff_data = mpi_diff.variables['difference_in_eto'][:]
mri_diff_data = mri_diff.variables['difference_in_eto'][:]
ukesm_diff_data = ukesm_diff.variables['difference_in_eto'][:]

mean_of_differences = (gfdl_diff_data + ipsl_diff_data + mpi_diff_data + mri_diff_data + ukesm_diff_data)/5
mean_of_differences_over_time = np.nanmean(mean_of_differences , axis=0)

lat = gfdl_diff.variables['lat'][:]
lon = gfdl_diff.variables['lon'][:]

# Create a plot
ax = plt.axes(projection=ccrs.PlateCarree())

plot = plt.contourf(lon, lat, mean_of_differences_over_time, cmap='rainbow', transform=ccrs.PlateCarree(), levels=np.linspace(0, 2.8, 20))
plt.colorbar(plot, ax=ax, shrink=0.8)

ax.set_extent([67.75, 97.75, 7.25, 37.25])
states.boundary.plot(ax=ax, linewidth=0.25, edgecolor='black')
plt.title('Mean of difference in Evapotranspiration (avg of 5GCMs) in 2021-2030')

plt.show()
