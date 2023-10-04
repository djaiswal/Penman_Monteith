from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature

# Open the NetCDF file
data = Dataset(r"D:\GFDL_ESM4_new\percentage_change_2021-2030",'r')

# Extract latitude, longitude, and time
lat = data.variables['latitude'][:]
lon = data.variables['longitude'][:]
time = data.variables['time'][:]

# Extract temperature data for all time steps
percent_change = data.variables['percentge_change'][:]

# Calculate the mean temperature over all time steps
mean_tas = np.mean(percent_change, axis=0)

# Create a plot
ax = plt.axes(projection=ccrs.PlateCarree())

plot = plt.contourf(lon, lat, mean_tas, cmap='viridis', transform=ccrs.PlateCarree())
plt.colorbar(plot, ax=ax, shrink=0.8)

ax.set_extent([67.75, 97.75, 7.25, 37.25])  # regional map (x0, x1, y0, y1)
ax.coastlines()
ax.gridlines(draw_labels=False)

ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.LAND, edgecolor='black')
#ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
#ax.add_feature(cartopy.feature.RIVERS)
plt.title('Percentage change in Evapotranspiration from 2021-2030 considering CO2')

plt.show()