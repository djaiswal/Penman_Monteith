# NOTE : This code calculates the Reference evapotranspiration for a decade (eg: 2021-2030, 2051-2060, or 2091-2100) for the ssp585 scenario.
#        The ETo values calculated without considering the effect of CO2 and ETo values calculated considering
#        the effect of CO2 are saved separately in a selected directory.

#        Enter the starting year of the decade for which reference evapotranspiration has to be calculated.
#        i.e, 2021 for the decade 2021-2030, 2051 for the decade 2051-2060, or 2091 for the decade 2091-2100.
#        Select the files containing the weather data from the first window. 
#        Select the parent folder where the ETo files for each GCM has to be saved. (For Eg: ETo_files)
#        Separate folders are created for different GCMs and respective files are saved there.
#        This creates a directory structure as follows:
# 
#           ETo_files
#                |--- GFDL-ESM4 (contains ETo files for GFDL-ESM4)
#                |--- IPSL-CM6A-LR (contains ETo files for IPSL-CM6A-LR)
#                |--- MPI-ESM1-2-HR (contains ETo files for MPI-ESM1-2-HR)
#                |--- MRI-ESM2-0 (contains ETo files for MRI-ESM2-0)
#                |--- UKESM1-0-LL (contains ETo files for UKESM1-0-LL)


import numpy as np
import netCDF4 as nc
import tkinter as tk
from tkinter import filedialog
import os
# import calculate_rn as calc_rn

# This part of the code creates the interface to browse the required files and directories.

file_paths = [""]*7
file_list = ['Relative humidity', 'Pressure',"Surface downwelling radiation",
             'Windspeed', 'Daily Max Temperature', 'Daily Min Temperature', 'Elevation']

# Get the initial conditions (starting year of decade under consideration)
start_year = input("Enter starting year: ")


def browse_file(file_index):
    file_path = filedialog.askopenfilename(
        filetypes=[(start_year+"files", "*"+start_year+"*"), ("Elevation files", "*elevation*"), ("NetCDF files", "*.nc") ]) # include file containing "elevation" in its name
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {file_list[file_index]} : {file_path}")

directory_path = ""

def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(
            text=f"Selected Directory to save ETo files : {directory_path}")

def finish():
    root.destroy()

root = tk.Tk()
root.title("File Browser")
root.geometry("800x600")  # Set the window size to 800x600

labels = []
buttons = []
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

root = tk.Tk()
root.title("Directory Browser to save ETo files")
root.geometry("800x200")  # Set the window size to 800x200
directory_button = tk.Button(
    root, text="Browse Directory to save ETo files", command=browse_directory)
directory_button.pack(pady=10)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=10)
finish_button = tk.Button(
    root, text="Finish selecting directory", command=finish)
finish_button.pack(pady=10)
root.mainloop()

# Identify the GCM from the file name (assuming that weather data files contain the name of GCM in their file name)
if 'gfdl' in file_paths[0]:
    gcm_name = 'GFDL-ESM4'
elif 'ipsl' in file_paths[0]:
    gcm_name = 'IPSL-CM6A-LR'
elif 'mpi' in file_paths[0]:
    gcm_name = 'MPI-ESM1-2-HR'
elif 'mri' in file_paths[0]:
    gcm_name = 'MRI-ESM2-0'
elif 'ukesm' in file_paths[0]:
    gcm_name = 'UKESM1-0-LL'

# Create a new directory to save the NetCDF files according to GCMs 7(if it is not present already)
# os.makedirs(output_parent_directory, exist_ok=True)
output_directory = os.path.join(directory_path, gcm_name)
os.makedirs(output_directory, exist_ok=True)

# This part of the code calculates net radiation from rsds,tasmax, tasmin, elevation and hurs and saves them as NetCDF files.

def calculate_rn(rsds, lat, tasmax, tasmin, ea, z, J):
    phi = (lat/180)*np.pi     # convert latitude in decimal degree to latitude in radian
    sigma = 4.903 * (10**-9)    # Steffan boltzmann constant in MJ K^-4 m^-2 day^-1
    # rs = rsds
    rns = 0.77*rsds     # net solar shortwave radiation 
    gsc = 0.0820    # solar constant in unit MJ m^-2 min^-1 
    dr = 1 + 0.033*np.cos(2*np.pi*J/365)    # inverse relative earth-sun distance
    delta = 0.409*np.sin((2*np.pi*J/365) - 1.39)
    omega_s = np.arccos((-np.tan(phi))*np.tan(delta))
    ra = (24*60/np.pi)*gsc*dr*(omega_s*np.sin(phi)*np.sin(delta) + np.cos(phi)*np.cos(delta)*np.sin(omega_s))
    # print(f"ra shape: {ra.shape}, rsds shape: {rsds.shape}")
    rso = (0.75 + 2*(10**-5)*z)*ra
    rnl = sigma*((tasmax**4 + tasmin**4)/2)*(0.34 -(0.14*((ea)**(1/2))))*(1.35*(rsds/rso) - 0.35)
    rn = rns - rnl
    return rn

def create_netrad_netcdf(file_name, output_directory, ref_dataset, start_year, rn_data):
    output_file = os.path.join(output_directory, file_name)
    dataset = nc.Dataset(output_file, 'w', format='NETCDF4')

    time_dim = len(ref_dataset.dimensions['time'])
    lat_dim = len(ref_dataset.dimensions['lat'])
    lon_dim = len(ref_dataset.dimensions['lon'])

    dataset.createDimension('time', time_dim)
    dataset.createDimension('lat', lat_dim)
    dataset.createDimension('lon', lon_dim)

    time_var = dataset.createVariable('time', np.float32, ('time',))
    lat_var = dataset.createVariable('lat', np.float32, ('lat', ))
    lon_var = dataset.createVariable('lon', np.float32, ('lon', ))

    if start_year == '2021':
        time_var.units = "days since 2021-1-1 00:00:00"
    elif start_year == '2091':
        time_var.units = "days since 2091-1-1 00:00:00"
    elif start_year == '2051':
        time_var.units = "days since 2051-1-1 00:00:00"
    lat_var.units = "degrees_north"
    lon_var.units = "degrees_east"

    radiation = dataset.createVariable(
        "net_radn", np.float32, ("time", "lat", "lon", )
    )

    radiation.units = "MJ "         # complete this line
    radiation.long_name = "Net radiation"

    time_var[:] = ref_dataset.variables["time"][:]
    lat_var[:] = ref_dataset.variables["lat"][:]
    lon_var[:] = ref_dataset.variables["lon"][:]
    radiation[:] = rn_data

    print(f"Net radiation data saved to {output_file}")
    dataset.close()
    return

# This part of the code computes the reference ETo and saves them as NetCDF files.

# Define constants for the Penman-Monteith equation
G = 0                                                                            # Soil heat flux
alpha = 0.23                                                                     # Albedo (typical value for grass reference crop)
gs_ref = 0.0061

# creating an array of size (3652, 61, 61) since we have 3652 timesteps, 61 latitude, 61 longitude
# co2_ref = np.full((3652, 61, 61), 330)

# Number of years including leap years
num_years = 10
days_per_year_2021 = [366 if year % 4 == 0 and (year % 100 != 0
                                           or year % 400 == 0)
                 else 365 for year in range(2021, 2031)]                           # Days per year, considering leap years

days_per_year_2091 = [366 if year % 4 == 0 and (year % 100 != 0
                                           or year % 400 == 0)
                 else 365 for year in range(2091, 2101)]

days_per_year_2051 = [366 if year % 4 == 0 and (year % 100 != 0
                                           or year % 400 == 0)
                 else 365 for year in range(2051, 2061)]                           # Days per year, considering leap years

# Function to calculate gs using the concentration of CO2
# def calculate_gs(start_year, gs_ref, co2_input):
def calculate_gs(start_year, gs_ref):

    # Creating arrays for conc of Co2
    if start_year == '2021':
        data1 = np.array([418.06, 421.33, 424.72, 428.22, 431.83,                   # CO2 concentrations in the years from 2021-2030
                          435.55, 439.38, 443.31, 447.36, 451.51])
        
        # Repeat the data for each day of the year
        data1_repeated_per_day = np.repeat(data1, days_per_year_2021)                    
        # Expand the data across all latitudes and longitudes
        data1_expanded = np.repeat(
            data1_repeated_per_day[:, np.newaxis, np.newaxis], 61 * 61, axis=-1)    
        # Reshape the data into the desired 3D shape
        data1_3dco2 = data1_expanded.reshape(3652, 61, 61)
        data_co2 = data1_3dco2

    elif start_year == '2091':
        data2 = np.array([1012.02, 1025.74, 1039.45, 1053.15, 1066.85,              # CO2 concentrations in the years from 2091-2100
                          1080.53, 1094.21, 1107.89, 1121.55, 1135.21])
        # Repeat the data for each day of the year
        data2_repeated_per_day = np.repeat(data2, days_per_year_2091)
        # Expand the data across all latitudes and longitudes
        data2_expanded = np.repeat(
            data2_repeated_per_day[:, np.newaxis, np.newaxis], 61 * 61, axis=-1)    
        # Reshape the data into the desired 3D shape
        data2_3dco2 = data2_expanded.reshape(3652, 61, 61)
        data_co2 = data2_3dco2
       
    elif start_year == '2051':
        # CO2 concentrations in the years from 2051-2060
        """
            2051  569.93
            2052  577.26
            2053  584.78
            2054  592.51
            2055  600.43
            2056  608.55
            2057  616.87
            2058  625.39
            2059  634.11
            2060  643.04
        """
        data3 = np.array([569.93, 577.26, 584.78, 592.51, 600.43,                   # CO2 concentrations in the years from 2051-2060
                          608.55, 616.87, 625.39, 634.11, 643.04])
        # Repeat the data for each day of the year
        data3_repeated_per_day = np.repeat(data3, days_per_year_2051)
        # Expand the data across all latitudes and longitudes
        data3_expanded = np.repeat(
            data3_repeated_per_day[:, np.newaxis, np.newaxis], 61 * 61, axis=-1)    
        # Reshape the data into the desired 3D shape
        data3_3dco2 = data3_expanded.reshape(3653, 61, 61)
        data_co2 = data3_3dco2

    # Calculate gs
    gs = gs_ref * (1 / (1 + 0.663 * ((data_co2 / 330) - 1)))
    return gs


# Function to calculate daily reference evapotranspiration (ET0)

def calculate_eto_with_CO2(tmean, wind_speed, relative_humidity, net_radiation, P, gs):

    # Calculate saturation vapor pressure (es) and actual vapor pressure (ea)
    es = 0.6108 * np.exp((17.27 * tmean) / (tmean + 237.3))
    ea = (relative_humidity / 100.0) * es

    # # Calculate net radiation
    # J = day % 365 + 1
    # Rn = calc_rn.calculate_rn(short_radiation, lat, temperature_max, temperature_min, ea, z, J)

    # Calculate the slope of the saturation vapor pressure curve (delta)
    delta = (4098 * es) / ((tmean + 237.3) ** 2)

    # Calculate the psychrometric constant (gamma)
    gamma = 0.000655*P

    # Calculate ET0 using the Penman-Monteith equation
    et0 = (((0.408 * delta * (net_radiation - G)) + (gamma * (900 / (tmean + 273)) * wind_speed * (es - ea))) /
           (delta + gamma * (1 + (0.0033 * (wind_speed/gs)))))

    return et0


def calculate_ETo_with_FAO(tmean, wind_speed, relative_humidity, net_radiation, P):

    es = 0.6108 * np.exp((17.27 * tmean) / (tmean + 237.3))
    ea = (relative_humidity / 100.0) * es

    # Calculate the slope of the saturation vapor pressure curve (delta)
    delta = (4098 * es) / ((tmean + 237.3) ** 2)

    # Calculate the psychrometric constant (gamma)
    gamma = 0.000655*P

    # Calculate ET0 using the Penman-Monteith equation
    et0 = (((0.408 * delta * (net_radiation - G)) + (gamma * (900 / (tmean + 273)) * wind_speed * (es - ea))) /
           (delta + gamma * (1 + (0.34 * (wind_speed)))))

    return et0    


# Function to create a new NetCDF file of ETo for individual GCMs and timeperiod

def create_netcdf(file_name, output_directory, time_dim, latitude_dim, longitude_dim, start_year, co2_input, ref_dataset, et0_values):
    # destination of output file
    output_file = os.path.join(output_directory, file_name)
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
    if start_year == '2021':
        time_var.units = 'days since 2021-1-1 00:00:00'
    elif start_year == '2091':
        time_var.units = 'days since 2091-1-1 00:00:00'
    elif start_year == '2051':
        time_var.units = 'days since 2051-1-1 00:00:00'
    latitude_var.units = 'degrees_north'
    longitude_var.units = 'degrees_east'
    # Define the variable for evapotranspiration
    evapotranspiration = dataset.createVariable(
        'evapotranspiration', np.float32, ('time', 'lat', 'lon',))

    # Define units and attributes for evapotranspiration
    evapotranspiration.units = 'mm/day'
    if co2_input == 'yes':
        evapotranspiration.long_name = 'Reference evapotranspiration considering the effect of CO2'
    else:
        evapotranspiration.long_name = 'Reference evapotranspiration with FAO-PM equation'

    # Assign data to variables
    time_var[:] = ref_dataset.variables['time'][:]
    latitude_var[:] = ref_dataset.variables['lat'][:]
    longitude_var[:] = ref_dataset.variables['lon'][:]
    evapotranspiration[:] = et0_values
    print(f"ET0 data saved to {output_file}")
    return

# Open the input NetCDF files
datahurs = nc.Dataset(file_paths[0], 'r')
datapressure = nc.Dataset(file_paths[1], 'r')
datarsds = nc.Dataset(file_paths[2], 'r')
datawindspeed = nc.Dataset(file_paths[3], 'r')
datatempmax = nc.Dataset(file_paths[4], 'r')
datatempmin = nc.Dataset(file_paths[5], 'r')
dataelevation = nc.Dataset(file_paths[6], 'r')

# Extract relevant variables from the input files

# in K in netcdf file, converted to degree celsius                              # Maximum temperature
temperature_max = datatempmax.variables['tasmax'][:] 
# in K in netcdf file , converted to degree celsius                             # Minimum temperature
temperature_min = datatempmin.variables["tasmin"][:]
tmean = ((temperature_min+temperature_max)/2) - 273
# surface wind speed in m/s                                                     # surface wind speed
wind_speed = datawindspeed.variables["sfcwind"][:]
# relative humidity in percentage                                               # relative humidity
relative_humidity = datahurs.variables["hurs"][:]
# in W/m^2 in netcdf file , converted to MJ/m^2                                 # net radiations
# net_radiation = datanetradn.variables["net_radn"][:] * (0.0864)
# in Pa in netcdf file , converted to KPa                                       # Atmospheric pressure
P = datapressure.variables["ps"][:] * (10**-3)
rsds_data = datarsds.variables["rsds"][:]
elevation_data = dataelevation.variables["Zenith_data"][:]
lat_data = datarsds.variables["lat"][:]
lat_data = lat_data[np.newaxis, :,  np.newaxis]
lat_data = np.repeat(lat_data, 61, axis=-1)  # Reshape lat_data to (1, 61, 61)
lat_data = np.repeat(lat_data, rsds_data.shape[0], axis=0)  # Reshape lat_data to (1, 61, 1)
lat_data = lat_data.reshape(rsds_data.shape)
days_per_year = [366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365 for year in range(int(start_year), int(start_year)+10)]
days = []
for ndays in days_per_year:
    for day in range(ndays):
        days.append((day+1))
days = np.array(days)
days_decade = np.repeat(days[:, np.newaxis, np.newaxis], 61*61, axis=-1)
days_decade = days_decade.reshape(rsds_data.shape)

def calculate_ea(temp_max, temp_min, rel_humidity):
    tmean = ((temp_min + temp_max)/2) - 273 # to convert into deg celcius
    es = 0.6108 * np.exp((17.27*tmean)/(tmean + 237.3))
    ea = (rel_humidity/100)*es
    return ea

ea_data = calculate_ea(temperature_max, temperature_min, relative_humidity)
rn_data = calculate_rn(rsds_data*0.0864, lat_data, temperature_max, temperature_min, ea_data, elevation_data, days_decade)
netrad_file = gcm_name + '_' + start_year + '-' + str(int(start_year)+9) + '_net_radiation.nc'
create_netrad_netcdf(netrad_file, output_directory, datahurs, start_year, rn_data)  # save the data as a NetCDF file

# Get dimensions from the input files
time_dim = len(datahurs.dimensions['time'])
latitude_dim = len(datahurs.dimensions['lat'])
longitude_dim = len(datahurs.dimensions['lon'])
ref_dataset = datahurs

file_name = gcm_name + '_eto_'
if start_year == '2021':
    file_name = file_name + '2021-2030_'
elif start_year == '2091':
    file_name = file_name + '2091-2100_'
elif start_year == '2051':
    file_name = file_name + '2051-2060_'

# Calculate ET0 for each day without considering the effect of CO2
co2_input = 'no'
et0_values_without_co2 = calculate_ETo_with_FAO(tmean,
                                       wind_speed, relative_humidity, rn_data, P)
# spatio_temporal_mean = np.mean(et0_values_without_co2, axis=(1, 2))  # Mean over lat and lon for each time step
overall_mean = np.mean(et0_values_without_co2, axis=(0, 1, 2))       # Mean over time, lat, and lon
print(f"Overall mean ET0 without CO2: {overall_mean:.2f} mm/day")  # Print the overall mean ET0
file_name_without_co2 = file_name + 'with_FAO-PM_eq.nc'
create_netcdf(file_name_without_co2, output_directory, time_dim, latitude_dim,                      # save the data as a NetCDF file
              longitude_dim, start_year, co2_input, ref_dataset, et0_values_without_co2)

# Calculate ET0 for each day without considering the effect of CO2
co2_input = 'yes'
gs_with_co2 = calculate_gs(start_year, gs_ref)
et0_values_with_co2 = calculate_eto_with_CO2(tmean,
                                    wind_speed, relative_humidity, rn_data, P, gs_with_co2)
overall_mean_with_co2 = np.mean(et0_values_with_co2, axis=(0, 1, 2))       # Mean over time, lat, and lon
print(f"Overall mean ET0 with CO2: {overall_mean_with_co2:.2f} mm/day")  # Print the overall mean ET0
file_name_with_co2 = file_name + 'with_CO2.nc'
create_netcdf(file_name_with_co2, output_directory, time_dim, latitude_dim,                        # save the data as a NetCDF file
              longitude_dim, start_year, co2_input, ref_dataset, et0_values_with_co2)

# Close the input and output files
datahurs.close()
datatempmax.close()
datatempmin.close()
datawindspeed.close()
datarsds.close()
dataelevation.close()
datapressure.close()
