import math
import pandas as pd                         # Read and nodify dataframe
import numpy as np
from netCDF4 import Dataset                 # Read netcdf files
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt             # to plot data
from datetime import datetime, timedelta    # manipulate date and time
import xarray as xr

# Open the first NetCDF file
ds_rlds = xr.open_dataset(r'D:\IITPKD\Main\5TH SEM\Tanushri Intern\WeatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rlds_ind_daily_2021_2030.nc',
                           engine = 'netcdf4')
print(ds_rlds)

# Open the second NetCDF file
ds_rsds = xr.open_dataset(r'D:\IITPKD\Main\5TH SEM\Tanushri Intern\WeatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rsds_ind_daily_2021_2030.nc',
                        engine = 'netcdf4')
print(ds_rsds.variables)
# Extract the hurs and ps variables from each dataset
rlds = ds_rlds['rlds']
rsds = ds_rsds['rsds']

# Add the variables together
net_radiation = (rlds + rsds)*(1 - 0.23)

# Update the fill values in the resulting net radiation variable
#fill_value = ds_rlds['rlds'].encoding['_FillValue']
# Set the fill value in the new variable
#net_radiation = net_radiation.where(~np.isnan(net_radiation), fill_value)


# Create a new dataset with the result
ds_net_radiation = xr.Dataset( 
    data_vars={
       'rn':xr.DataArray(net_radiation,  
                          coords = {'time' : ds_rlds['time'], 'lat' : ds_rlds['lat'], 'lon' : ds_rlds['lon']},
                          attrs = ds_rlds['rlds'].attrs)})

print(ds_net_radiation.variables)
print(ds_rlds.variables)

# Save the new dataset to a NetCDF file
ds_net_radiation.to_netcdf('net_radiation.nc')
