import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
from netCDF4 import Dataset
import numpy as np

shapefile_path = r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\india_state_boundary_projected.shp"
states = gpd.read_file(shapefile_path)
states = states.to_crs(ccrs.PlateCarree().proj4_init)

gfdl_diff = Dataset(r"E:\ISIMIP Climate Data\eto_masked_changed\difference_netcdf_files\masked_gfdl_et0_2021-2030_difference.nc",'r')
ipsl_diff = Dataset(r"E:\ISIMIP Climate Data\eto_masked_changed\difference_netcdf_files\masked_ipsl_et0_2021-2030_difference.nc",'r')
mpi_diff = Dataset(r"E:\ISIMIP Climate Data\eto_masked_changed\difference_netcdf_files\masked_mpi_et0_2021-2030_difference.nc",'r')
mri_diff = Dataset(r"E:\ISIMIP Climate Data\eto_masked_changed\difference_netcdf_files\masked_mri_et0_2021-2030_difference.nc",'r')
ukesm_diff = Dataset(r"E:\ISIMIP Climate Data\eto_masked_changed\difference_netcdf_files\masked_ukesm_et0_2021-2030_difference.nc",'r')

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

plot = plt.contourf(lon, lat, mean_of_differences_over_time, cmap='YlGn', transform=ccrs.PlateCarree())
plt.colorbar(plot, ax=ax, shrink=0.8)

ax.set_extent([67.75, 97.75, 7.25, 37.25])
states.boundary.plot(ax=ax, linewidth=0.25, edgecolor='black')
plt.title('Mean of difference in Evapotranspiration (avg of 5GCMs) in 2021-2030')

plt.show()
