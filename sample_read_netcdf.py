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
ps = Dataset(r'D:\IITPKD\Main\5TH SEM\Tanushri Intern\WeatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_ps_ind_daily_2021_2030.nc')

# Printing the variables in netcdf file
print(hurs)
print(hurs.variables.keys())
# Reading details of latitude variable
lats = hurs.variables['lat']
print(lats)
# Assign variable to latitude and read all the values with [:]
lats = hurs.variables['lat'][:]
print(lats)
# Read first 10 values of latitude
lats = hurs.variables['lat'][:10]
print(lats)
# Reading longitude and printing values
lons = hurs.variables['lon'][:]
print(lons)
# Reading time and printing values
time = hurs.variables['time'][:]
print(time)
# to know th units of time
unit_time = hurs.variables['time'].units
print(unit_time)
# to know th units of rh
unitrh = hurs.variables['hurs'].units
print(unitrh)
