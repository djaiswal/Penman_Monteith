import numpy as np
import netCDF4 as nc

# Define constants for the Penman-Monteith equation

G = 0                               # Soil heat flux 
alpha = 0.23                        # Albedo (typical value for grass reference crop) 
gs_ref=0.0061
co2_ref=np.full((3652, 61, 61), 330)   #creating an array of same size as netcdf

#get the initial conditions
co2_input=input("consider co2 concentration?")
start_year=input("enter starting year: ")

#creating arrays for conc of Co2
num_years = 10                                                        # Number of years including leap years
data1 = np.array([418.06,421.33,424.72,428.22,431.83,435.55,439.38,443.31,447.36,451.51])
data2 = np.array([1012.02,1025.74,1039.45,1053.15,1066.85,1080.53,1094.21,1107.89,1121.55,1135.21])

days_per_year = [366 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 365 for year in range(2021, 2031)]          # Days per year, considering leap years

# Repeat the data for each day of the year
data1_repeated_per_day = np.repeat(data1, days_per_year)
data2_repeated_per_day = np.repeat(data2, days_per_year)
# Expand the data across all latitudes and longitudes
data1_expanded = np.repeat(data1_repeated_per_day[:, np.newaxis, np.newaxis], 61 * 61, axis=-1)
data2_expanded = np.repeat(data2_repeated_per_day[:, np.newaxis, np.newaxis], 61 * 61, axis=-1)
# Reshape the data into the desired 3D shape
data1_3dco2 = data1_expanded.reshape(3652, 61, 61)
data2_3dco2 = data2_expanded.reshape(3652, 61, 61)

#select whether to include co2 or not
if co2_input=="yes" and start_year=='2021':
    data_co2=data1_3dco2
elif co2_input=="yes" and start_year=='2091':
    data_co2=data2_3dco2
else:
    data_co2=co2_ref

gs = gs_ref * (1 / (1 + 0.663 * ((data_co2/ co2_ref) - 1))) 

# Function to calculate daily reference evapotranspiration (ET0) 

def calculate_et0(temperature_max, temperature_min, wind_speed, relative_humidity, long_radiation,short_radiation,P):

    tmean = ((temperature_max + temperature_min) / 2.0)                            # Calculate daily mean temperature (Tmean) in Celsius
    Rn=(long_radiation + short_radiation)*(0.77)
    
    es = 0.6108 * np.exp((17.27 * tmean) / (tmean + 237.3))        # Calculate saturation vapor pressure (es) and actual vapor pressure (ea)
    ea = (relative_humidity / 100.0) * es
           
    delta = (4098 * es) / ((tmean + 237.3) ** 2)                   # Calculate the slope of the saturation vapor pressure curve (delta)

    gamma = 0.000655*P                       # Calculate the psychrometric constant (gamma)

    et0 = ((0.408 * delta * (Rn - G) + (gamma * (900 / (tmean + 273)) * wind_speed * (es - ea))) /         # Calculate ET0 using the Penman-Monteith equation
           (delta + gamma * (1 + 0.035 * wind_speed)/gs))

    return et0

    


if start_year=='2021':
# Open the input NetCDF files 

    datahurs = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_hurs_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

    datatempmax = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasmax_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

    datatempmin = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasmin_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

    datawindspeed = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_sfcwind_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

    datalongradn = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rlds_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

    datashortradn = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rsds_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

    datapressure = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_ps_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc",'r')

elif start_year=='2091':

    datahurs = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_hurs_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')

    datatempmax = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasmax_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')

    datatempmin = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasmin_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')

    datawindspeed = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_sfcwind_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')

    datalongradn = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rlds_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')

    datashortradn = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rsds_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')

    datapressure = nc.Dataset(r"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_ps_lat7.25to37.25lon67.75to97.75_daily_2091_2100.nc",'r')


# Extract relevant variables from the input files
#"E:\ISIMIP Climate Data\GFDL-ESM4\7.25N 37.25N 67.75E 97.75E\gfdl-esm4_r1i1p1f1_w5e5_ssp585_hurs_lat7.25to37.25lon67.75to97.75_daily_2021_2030.nc"

temperature_max = datatempmax.variables['tasmax'][:] - 273  #in K in netcdf file , converted to degree celsius     #Maximum temperature
temperature_min = datatempmin.variables["tasmin"][:] - 273  #in K in netcdf file , converted to degree celsius     #Minimum temperature
wind_speed = datawindspeed.variables["sfcwind"][:]                                                                 #surface wind speed
relative_humidity = datahurs.variables["hurs"][:]                                                                  #relative humidity
long_radiation = datalongradn.variables["rlds"][:]* (0.0864)      #in W/m^2 in netcdf file , converted to MJ/m^2       #longwave radiation
short_radiation = datashortradn.variables["rsds"][:]* (0.0864)    #in W/m^2 in netcdf file , converted to MJ/m^2       #shortwave radiations
P = datapressure.variables["ps"][:]* (10**-3)                     #in Pa in netcdf file , converted to KPa             #Atmospheric pressure 
latitude = datahurs.variables["lat"][:]                                                                            #latitudes
longitude = datahurs.variables["lon"][:]                                                                           #longitudes

# Calculate ET0 for each day
et0_values = calculate_et0(temperature_max, temperature_min, wind_speed, relative_humidity, long_radiation, short_radiation,P)

# Get dimensions from the input files
time_dim = len(datahurs.dimensions['time'])
latitude_dim = len(datahurs.dimensions['lat'])
longitude_dim = len(datahurs.dimensions['lon'])

# Create a new NetCDF file
output_file = r"E:\ISIMIP Climate Data\eto_equation_changed\gfdl_esm4_changed_et0_2021-2030_with_CO2.nc"
dataset = nc.Dataset(output_file, 'w', format='NETCDF4')

# Define dimensions in the NetCDF file
time = dataset.createDimension('time', time_dim)
latitude = dataset.createDimension('lat', latitude_dim)
longitude = dataset.createDimension('lon', longitude_dim)

# Create variables for time, latitude, and longitude
time_var = dataset.createVariable('time', np.float32, ('time',))
latitude_var = dataset.createVariable('lat', np.float32, ('lat',))
longitude_var = dataset.createVariable('lon', np.float32, ('lon',))

# Define units for variables
time_var.units = 'days since 2021-1-1 00:00:00'
latitude_var.units = 'degrees_north'
longitude_var.units = 'degrees_east'
# Define the variable for evapotranspiration
evapotranspiration = dataset.createVariable('evapotranspiration', np.float32, ('time', 'lat', 'lon',))

# Define units and attributes for evapotranspiration
evapotranspiration.units = 'mm/day'
evapotranspiration.long_name = 'Reference evapotranspiration considering the effect of CO2'
# without
# Assign data to variables
time_var[:] = np.arange(time_dim)
latitude_var[:] = datahurs.variables['lat'][:]
longitude_var[:] = datahurs.variables['lon'][:]
evapotranspiration[:] = et0_values


# Close the input and output files
datahurs.close()
datatempmax.close()
datatempmin.close()
datawindspeed.close()
datalongradn.close()
datashortradn.close()
datapressure.close()


print(f"ET0 data saved to {output_file}")
