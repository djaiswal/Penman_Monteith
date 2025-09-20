# NOTE : This code calculates the daily net radiation for a decade (2021-2030, 2051-2060, or 2091-2100). 

#        Select the files containing the weather data from the first window. 
#        Select the parent folder where the files with data of net radiation has to be saved(For Eg: ETo_files)

import netCDF4 as nc
import tkinter as tk
from tkinter import filedialog
import os
import numpy as np

filepaths = [""]*5
file_list = ["Surface downwelling radiation", "Maximum temperatute", 
             "Minimum temperature", "Relative humidity", "Elevation"]

def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        filepaths[file_index] = file_path
        labels[file_index].config(
            text=f"{file_list[file_index]} : {file_path}"
        )

def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(text = f"Selected directory : {directory_path}")

def finish():
    root.destroy()
root = tk.Tk()
root.title("File browser")
root.geometry("800x600")

labels = []
buttons = []
for i in range(len(filepaths)):
    button = tk.Button(
        root, text=f"Browse file for {file_list[i]}", command = lambda i=i:browse_file(i))
    button.pack(pady = 5)
    buttons.append(button)

    label = tk.Label(root, text=f"{file_list[i]} : ")
    label.pack(pady = 5)
    labels.append(label)

finish_button = tk.Button(root, text = "Finish selecting files", command = finish)
finish_button.pack(pady = 5)
root.mainloop()

root = tk.Tk()
root.title("Browser")
root.geometry("800x200")
directory_button = tk.Button(
    root, text="Browse directory to save net radiation files", command=browse_directory
)

directory_button.pack(pady = 10)
directory_label = tk.Label(root, text= "Selected Directory: ")
directory_label.pack(pady = 10)
finish_button = tk.Button(
    root, text="Finish selecting directory", command = finish
)
finish_button.pack(pady = 10)
root.mainloop()

if 'gfdl' in filepaths[0]:
    gcm_name = 'GFDL-ESM-04'
elif 'ipsl' in filepaths[0]:
    gcm_name = "IPSL-CM6A-LR"
elif 'mpi' in filepaths[0]:
    gcm_name = "MPI-ESM1-2-HR"
elif 'mri' in filepaths[0]:
    gcm_name = "MRI-ESM2-0"
elif 'ukesm' in filepaths[0]:
    gcm_name = 'UKESM1-0-LL'

#Create a new netcdf file to save the net radiation
global num_days_in_decade
num_days_in_decade = 3652 
if '2021' in filepaths[0]:
    start_year = 2021
    num_days_in_decade = 3652
elif '2091' in filepaths[0]:
    start_year = 2091
    num_days_in_decade = 3653
elif '2051' in filepaths[0]:
    start_year = 2051
    num_days_in_decade = 3652

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
    rso = (0.75 + 2*(10**-5)*z)*ra
    rnl = sigma*((tasmax**4 + tasmin**4)/2)*(0.34 -(0.14*((ea)**(1/2))))*(1.35*(rsds/rso) - 0.35)
    rn = rns - rnl
    return rn

def create_netcdf(file_name, output_directory, ref_dataset, start_year, rn_data):
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

    if start_year == 2021:
        time_var.units = "days since 2021-1-1 00:00:00"
    elif start_year == 2091:
        time_var.units = "days since 2091-1-1 00:00:00"
    elif start_year == 2051:
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


datarsds = nc.Dataset(filepaths[0], 'r')
datatempmax = nc.Dataset(filepaths[1], 'r')
datatempmin = nc.Dataset(filepaths[2], 'r')
datahurs = nc.Dataset(filepaths[3], 'r')
dataelevation = nc.Dataset(filepaths[4], 'r')

rsds_data = datarsds.variables['rsds'][:]
temp_max = datatempmax.variables['tasmax'][:]
temp_min = datatempmin.variables['tasmin'][:]
rel_humidity = datahurs.variables['hurs'][:]
elev_data = dataelevation.variables['Zenith_data'][:]

lat_data = datarsds.variables["lat"][:]
lat_data = lat_data[np.newaxis, :, np.newaxis]
lat_data = np.repeat(lat_data, 61, axis=-1)
lat_data = np.repeat(lat_data, num_days_in_decade, axis=0)
lat_data = lat_data.reshape(num_days_in_decade, 61, 61)


start_decade = start_year
end_decade = start_year + 9
days_per_year = []
for year in range(start_decade, end_decade + 1):
    # 2100 is not a leap year, all other years divisible by 4 are leap years
    if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) and year != 2100:
        days_per_year.append(366)
    else:
        days_per_year.append(365)
days = []
for ndays in days_per_year:
    for day in range(ndays):
        days.append((day+1))
days = np.array(days)
days_decade = np.repeat(days[:, np.newaxis, np.newaxis], 61*61, axis=-1)
days_decade = days_decade.reshape(num_days_in_decade, 61, 61)

def calculate_ea(temp_max, temp_min, rel_humidity):
    tmean = ((temp_min + temp_max)/2) - 273 # to convert into deg celcius
    es = 0.6108 * np.exp((17.27*tmean)/(tmean + 237.3))
    ea = (rel_humidity/100)*es
    return ea

ea_data = calculate_ea(temp_max, temp_min, rel_humidity)
rn_data = calculate_rn(rsds_data, lat_data, temp_max, temp_min, ea_data, elev_data, days_decade)

file_name = gcm_name + '_' + str(start_year) + '_' + str(start_year+9)+'_net_radiation.nc'
create_netcdf(file_name, directory_path, datarsds, start_year, rn_data)

datarsds.close()
datatempmax.close()
datatempmin.close()
datahurs.close()
dataelevation.close()


