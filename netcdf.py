import netCDF4 as nc
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

data = Dataset(r'C:\Python\Project_1\MERRA2_100.inst1_2d_lfo_Nx.19800104.SUB.nc')
print(data)
data.variables.keys()
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
time = data.variables['time'][:]
press = data.variables['PS'][:]
temp = data.variables['TLML'][:]
rh = data.variables['QLML'][:]
# making a dummy formula here
dum_formula = (press + temp) * rh
print(lons)
mp = Basemap(projection='merc',
             llcrnrlon=-50.0,
             llcrnrlat=10.0,
             urcrnrlon=-95.0,
             urcrnrlat=45.0,
             resolution='i')

lon, lat = np.meshgrid(lons, lats)
x, y = mp(lon, lat)
c_scheme = mp.pcolor(x, y, np.squeeze(dum_formula[0, :, :]), cmap='jet')
mp.drawcoastlines()
mp.drawstates()
mp.drawcounties()

parallels = np.arange(10, 50, 5.)  # define latitude lines every 5 degrees from 10N to 50N
meridians = np.arange(-95, -50, 5.)  # define longitude lines every 5 degrees from 95W to 50W
mp.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)

for meridian in meridians:
    lon, lat = mp(meridian, 10)  # Adjust latitude value for placing the label
    plt.annotate(str(meridian), xy=(lon, lat), xytext=(5, -5), textcoords='offset points', fontsize=8,
                 ha='center', va='top', rotation='horizontal')

plt.title('Dummy Formula')
plt.grid(True)
plt.show()

new_dataset = nc.Dataset('new_file.nc', 'w')  # Create a new NetCDF file in write mode

# Create dimensions for latitude, longitude, and time
new_dataset.createDimension('lat', len(lats))
new_dataset.createDimension('lon', len(lons))
new_dataset.createDimension('time', None)  # None allows unlimited dimension size for time

# Create variables for latitude, longitude, time, and the dummy formula
lat_var = new_dataset.createVariable('lat', 'f4', ('lat',))
lon_var = new_dataset.createVariable('lon', 'f4', ('lon',))
time_var = new_dataset.createVariable('time', 'f8', ('time',))
formula_var = new_dataset.createVariable('dummy_formula', 'f4', ('time', 'lat', 'lon',))

# Assign attributes to the variables
lat_var.units = 'degrees_north'
lon_var.units = 'degrees_east'
time_var.units = 'hours since 1970-01-01 00:00:00'
formula_var.units = 'dummy units'

# Assign the latitude, longitude, and time data to the corresponding variables
lat_var[:] = lats
lon_var[:] = lons
time_var[:] = time

formula_data = (press + temp) * rh

formula_var[:] = formula_data

new_dataset.close()


data.close()


import netcdf as nc 
from netCDF4 import Dataset 

data = Dataset(r'C:\Python\new_file.nc')
print(data)