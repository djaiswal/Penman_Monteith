import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
from netCDF4 import Dataset
import numpy as np
import os
import tkinter as tk
from tkinter import filedialog
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

gcm_name = input("Enter GCM name : ")
# Load shapefile with boundaries of Indian states
os.environ['SHAPE_RESTORE_SHX'] = 'YES'
shapefile_path = ""

file_paths = [""]*2
file_list = ['2021-2030', '2091-2100']
def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {file_list[file_index]} : {file_path}")

def browse_shapefile():
        global shapefile_path
        shapefile_path = filedialog.askopenfilename(filetypes=[("Shapefiles", "*.shp"), ("All files", "*.*")])
        label_1.config(text=f"Path to shapefile : {shapefile_path}")


def finish():
    root.destroy()

root = tk.Tk()
root.title("File Browser")
root.geometry("800x300")  # Set the window size to 800x1000

labels = []
buttons = []
button_1 = tk.Button(root, text="Browse shapefile", command=browse_shapefile)
button_1.pack(pady=5)
label_1 = tk.Label(root, text="Shapefile : ")
label_1.pack(pady=10)

for i in range(len(file_paths)):

    button = tk.Button(
        root, text=f"Browse File for {file_list[i]}", command=lambda i=i: browse_file(i))
    button.pack(pady=5)
    buttons.append(button)

    label = tk.Label(root, text=f"{file_list[i]}: ")
    label.pack(pady=5)
    labels.append(label)

finish_button = tk.Button(root, text="Finish Selecting Files", command=finish)
finish_button.pack(pady=5)
root.mainloop()

mean_et = []

for file_path in file_paths:
    data = Dataset(file_path, 'r')
    et = data.variables['difference_in_eto'][:]
    mean_et.append(np.mean(et, axis=0))
states = gpd.read_file(shapefile_path)
states = states.set_crs('EPSG:4326')
states = states.to_crs(ccrs.PlateCarree().proj4_init)

max_value = max(np.nanmax(mean_et[0]), np.nanmax(mean_et[1]))
# print(max_value)

data1 = Dataset(file_paths[0], 'r')
lat = data1.variables['lat'][:]
lon = data1.variables['lon'][:]
time = data1.variables['time'][:]

original_cmap = plt.colormaps.get_cmap('tab20b')

def get_partial_colormap(cmap, start=0.0, end=1.0, n=256):
    cmap_array = cmap(np.linspace(start, end, n))
    new_cmap = LinearSegmentedColormap.from_list('partial_cmap', cmap_array)
    return new_cmap

def interpolate_cmap(cmap, num_colors=256):
    # Get the list of original colors
    colors = cmap(np.linspace(0, 1, cmap.N))
    
    # Create an interpolated version of the colormap
    interpolated_colors = np.zeros((num_colors, 4))
    for i in range(4):  # Iterate over RGBA
        interpolated_colors[:, i] = np.interp(np.linspace(0, 1, num_colors), np.linspace(0, 1, cmap.N), colors[:, i])
    
    # Create a new colormap
    new_cmap = ListedColormap(interpolated_colors)
    return new_cmap

# Increase the number of color variants
num_intervals = 44
# num_intervals = int(num_intervals)
new_cmap = interpolate_cmap(original_cmap, num_colors=num_intervals)
new_cmap = get_partial_colormap(new_cmap, start=0.075, end=1.0)

fig, axes = plt.subplots(1, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(20, 10))
levels = np.linspace(0, 3.000, num_intervals)
for i, ax in enumerate(axes.flat):
    plot = ax.contourf(lon, lat, mean_et[i], cmap=new_cmap, transform=ccrs.PlateCarree(), levels=levels)
    ax.set_extent([67.75, 97.75, 7.25, 37.25])  
    states.boundary.plot(ax=ax, linewidth=0.25, edgecolor='black')
    gl1 = ax.gridlines(draw_labels=True, alpha=0.25, linestyle='--')
    gl1.top_labels = False
    gl1.right_labels = False
    titles = ['2021-2030', '2091-2100']
    ax.set_title(titles[i])

# Create colorbar0.735, 0.25, 0.02, 0.5
# cax0 = fig.add_axes([0.735, 0.25, 0.02, 0.5])  # Position of colorbar for subplot 1
cbar0 = plt.colorbar(plot, ax = axes, orientation='horizontal', label='Reference Evapotranspiration (mm/day)', shrink=0.8, pad=0.05, aspect=20, ticks=levels)
cbar0.ax.tick_params(rotation=90)

plt.suptitle(gcm_name.format(max_value), fontsize=18)
plt.subplots_adjust()
plt.savefig('mean_differences_plot.png', dpi=300, bbox_inches='tight')
plt.show()
