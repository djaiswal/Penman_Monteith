from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature

# Define the paths to your four NetCDF files to be plotted
file_paths = [
    r"D:\GFDL_ESM4_new\evapotranspiration_with_co2_2091-2100.nc",
    r"D:\GFDL_ESM4_new\evapotranspiration_without_co2_2091-2100.nc",
    r"D:\GFDL_ESM4_new\evapotranspiration_with_co2_2021-2030.nc",
    r"D:\GFDL_ESM4_new\evapotranspiration_without_co2_2021-2030.nc",
]

# Create an empty list to store the mean data for each dataset
mean_et = []

# Loop through the NetCDF files to extract and calculate mean data
for file_path in file_paths:
    data = Dataset(file_path, 'r')
    et = data.variables['evapotranspiration'][:]
    mean_et.append(np.mean(et, axis=0))

# Find the maximum value in all datasets to normalize them
max_value = max(np.max(mean_et[0]), np.max(mean_et[1]), np.max(mean_et[2]), np.max(mean_et[3]))
data1 = Dataset(r"D:\GFDL_ESM4_new\evapotranspiration_with_co2_2091-2100.nc", 'r')
lat = data1.variables['latitude'][:]
lon = data1.variables['longitude'][:]
time = data1.variables['time'][:]

# Create a 2x2 grid of subplots
fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(16, 10))

# Create subplots and add data
for i, ax in enumerate(axes.flat):
    plot = ax.contourf(lon, lat, mean_et[i], cmap='RdYlBu_r', transform=ccrs.PlateCarree(), levels=np.linspace(0, max_value, 11))
    ax.set_extent([67.75, 97.75, 7.25, 37.25])  # regional map (x0, x1, y0, y1)
    ax.coastlines()
    ax.gridlines(draw_labels=False)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND, edgecolor='black')
    titles = ['Evapotranspiration considering CO2 2091-2100','Evapotranspiration without considering CO2 2091-2100', 'Evapotranspiration considering CO2 2021-2030', 'Evapotranspiration without considering CO2 2021-2030']
    ax.set_title(titles[i])

# Create colorbars for each subplot
cax0 = fig.add_axes([0.92, 0.55, 0.02, 0.3])  # Position of colorbar for subplot 1
cax1 = fig.add_axes([0.92, 0.15, 0.02, 0.3])  # Position of colorbar for subplot 2
cax2 = fig.add_axes([0.47, 0.55, 0.02, 0.3])  # Position of colorbar for subplot 3
cax3 = fig.add_axes([0.47, 0.15, 0.02, 0.3])  # Position of colorbar for subplot 4

cbar0 = plt.colorbar(plot, cax=cax0, orientation='vertical', label='Evapotranspiration (in mm/day)', extend='both')
cbar1 = plt.colorbar(plot, cax=cax1, orientation='vertical', label='Evapotranspiration (in mm/day)', extend='both')
cbar2 = plt.colorbar(plot, cax=cax2, orientation='vertical', label='Evapotranspiration (in mm/day)', extend='both')
cbar3 = plt.colorbar(plot, cax=cax3, orientation='vertical', label='Evapotranspiration (in mm/day)', extend='both')

plt.suptitle('Comparison of Mean Evapotranspiration'.format(max_value),fontsize=20)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
