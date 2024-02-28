import math
import pandas as pd                         # Read and nodify dataframe
import numpy as np
from netCDF4 import Dataset                 # Read netcdf files
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt             # to plot data
from datetime import datetime, timedelta    # manipulate date and time
#import nctoolkit as nc

# Reading netcdf files, 'r' stands for reading mode
hurs = Dataset(r'D:\IITPKD\Main\5TH SEM\Tanushri Intern\WeatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_hurs_ind_daily_2021_2030.nc')
lats = hurs.variables['lat'][:]
lons = hurs.variables['lon'][:]
time = hurs.variables['time'][:]

# Setting Boundary
latbound = [6.75,37.75]
lonbound = [67.75,97.75]

# To get the index of the 4 points of the boundary
lat_lb = np.argmin(abs(lats - latbound[0]))
lat_ub = np.argmin(abs(lats - latbound[1]))
lon_lb = np.argmin(abs(lons - lonbound[0]))
lon_ub = np.argmin(abs(lons - lonbound[1]))

# To get the lat lon of the subset file
lats_sub = hurs.variables['lat'][lat_ub:lat_lb]
lons_sub = hurs.variables['lon'][lon_lb:lon_ub]
print(lats_sub)
print(lons_sub)
