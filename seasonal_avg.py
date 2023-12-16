
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, timedelta
import geopandas as gpd
import os
import glob

# Directory containing NetCDF files
directory_path = r'E:\ISIMIP Climate Data\eto_masked_changed\eto_files\ukesm_eto'

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
    base_date = datetime.strptime('2021-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
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

    print("mean_timeandspace_winter_eto", mean_timeandspace_winter_eto)
    print("mean_timeandspace_pre_monsoon_eto", mean_timeandspace_pre_monsoon_eto)
    print("mean_timeandspace_monsoon_eto",mean_timeandspace_monsoon_eto)
    print("mean_timeandspace_post_monsoon_eto", mean_timeandspace_post_monsoon_eto )
        
    mean_et=[]
    mean_et.append(mean_winter_eto)
    mean_et.append(mean_pre_monsoon_eto)
    mean_et.append(mean_monsoon_eto)
    mean_et.append(mean_post_monsoon_eto)

    lat = dataset.variables['lat'][:]
    lon = dataset.variables['lon'][:]
    time = dataset.variables['time'][:]

    shapefile_path = r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\india_state_boundary_projected.shp"
    states = gpd.read_file(shapefile_path)
    states = states.to_crs(ccrs.PlateCarree().proj4_init)

    max_value =max(np.nanmax(mean_et[0]), np.nanmax(mean_et[1]), np.nanmax(mean_et[2]), np.nanmax(mean_et[3]))
    fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(16, 10))
    plt.subplots_adjust(left=0.275, right=0.725, top=0.9, bottom=0.075, wspace=0.02, hspace=0.2)
    for i, ax in enumerate(axes.flat):
        plot = ax.contourf(lon, lat, mean_et[i], cmap='pink_r', transform=ccrs.PlateCarree(), levels=np.linspace(0, 9.5, 20))
        ax.set_extent([67.75, 97.75, 7.25, 37.25])  
        states.boundary.plot(ax=ax, linewidth=0.25, edgecolor='black')
        titles = ['mean_winter_eto','mean_pre_monsoon_eto', 'mean_monsoon_eto', 'mean_post_monsoon_eto']
        ax.set_title(titles[i])

    # Create colorbar
    cax0 = fig.add_axes([0.735, 0.25, 0.02, 0.5])  # Position of colorbar for subplot 1
    cbar0 = plt.colorbar(plot, cax=cax0, orientation='vertical', label='Evapotranspiration (in mm/day)', extend='both')

    plt.suptitle(file_path[file_path.find('1_')+2:file_path.find('.nc')].format(max_value),fontsize=18)
    plt.show()
    dataset.close()