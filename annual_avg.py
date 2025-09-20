
# NOTE : This code calculates annual avg of ETo from daily ETo for both the timeperiods and five GCMs.
#        The calculated values are stored in an excel sheet. 

#        Select the directory containing NetCDF files of daily Eto from selection window.
#        Suppose the directory structure is as follows:
#           ETo_files
#                |--- GFDL-ESM4
#                |--- IPSL-CM6A-LR
#                |--- MPI-ESM1-2-HR
#                |--- MRI-ESM2-04-0
#                |--- UKESM1-0-LL
#
#        Select the directory 'ETo_files' in the selection window.
#
#        The second input required is the directory to which
#        annual average of Eto has to be saved.


import netCDF4 as nc
import numpy as np
from datetime import datetime
import pandas as pd
import os
import glob
import tkinter as tk
from tkinter import filedialog
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# This part of the code creates the interface to select the required directories
def browse_directory():
    global directory_path_var_1
    path = filedialog.askdirectory()
    if path:
        directory_path_var_1.set(path)
        directory_label.config(text=f"Selected Directory with ETo files : {directory_path_var_1.get()}")

def browse_directory_1():
    global directory_path_var_2
    path2 = filedialog.askdirectory()
    if path2:
        directory_path_var_2.set(path2)
        directory_label_1.config(text=f"Selected Directory with ETo files : {directory_path_var_2.get()}")
def finish():
    root.destroy()

# Create the main window
root = tk.Tk()
root.title("Directory Browser ")
root.geometry("800x300")  # Set the window size to 800x200
directory_path_var_1 = tk.StringVar()
directory_path_var_2 = tk.StringVar()
directory_button = tk.Button(root, text="Browse Directory with daily data of ETo", command=browse_directory)
directory_button.pack(pady=10)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=10)
directory_button_1 = tk.Button(root, text="Browse Directory to save annual avg of ETo", command=browse_directory_1)
directory_button_1.pack(pady=10)
directory_label_1 = tk.Label(root, text="Selected Directory: ")
directory_label_1.pack(pady=10)
finish_button = tk.Button(root, text="Finish selecting directory", command=finish)
finish_button.pack(pady=10)
root.mainloop()

directory_path_var_1 = directory_path_var_1.get()
directory_path_var_2 = directory_path_var_2.get()

# This part of the code calculates the annual average of reference evapotranspiration
# Get the list of files and directories in the specified path
files_list = os.listdir(directory_path_var_1)
print(f"Files in the directory: {files_list}")
# Define the pattern to match NetCDF files
file_pattern = '*.nc'

# Create a list of file paths for all Netcdf files in files_list
file_paths_dir = [os.path.join(directory_path_var_1, f) for f in files_list]
file_paths = []
for directory_path in file_paths_dir:
    # Use glob to find all files matching the pattern
    if os.path.isdir(directory_path):
        file_paths.extend(glob.glob(os.path.join(directory_path, file_pattern)))
    else:
        print(f"Skipping {directory_path}: Not a directory.")
        continue

# file_paths = glob.glob(os.path.join(directory_path_var_1, file_pattern))

# Initialize an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['GCM','CO2', 'Year', 'Mean Ref. ETo'])

start_year = '20'
# Iterate through each file
for file_path in file_paths:
    if ('gfdl' in file_path) or ('GFDL' in file_path):
        gcm_name = 'GFDL-ESM4'
    elif ('ipsl' in file_path) or ('IPSL' in file_path):
        gcm_name = 'IPSL-CM6A-LR'
    elif ('mpi' in file_path) or ('MPI' in file_path):
        gcm_name = 'MPI-ESM1-2-HR'
    elif ('mri' in file_path) or ('MRI' in file_path):
        gcm_name = 'MRI-ESM2-0'
    elif ('ukesm' in file_path) or ('UKESM' in file_path):
        gcm_name = 'UKESM1-0-LL'
    dataset = nc.Dataset(file_path, 'r')
    if 'evapotranspiration' not in dataset.variables:
        print(f"Skipping {file_path}: 'evapotranspiration' variable not found.")
        dataset.close()
        continue

    time_var = dataset.variables['time']
    eto_var = dataset.variables['evapotranspiration'][:]

    # Assuming time is in a standard datetime format, you can use the following to extract the year
    if '2021' in file_path:
        start_year = "2021"
        base_date = "2021-01-01 00:00:00"
    elif '2091' in file_path:
        start_year = "2091"
        base_date = "2091-01-01 00:00:00"
    elif '2051' in file_path:
        start_year = "2051"
        base_date = "2051-01-01 00:00:00"
    base_datetime = nc.num2date(0, units=time_var.units, calendar=getattr(time_var, 'calendar', 'gregorian'))
    days_since_base_date = (base_datetime - datetime.strptime(base_date, "%Y-%m-%d %H:%M:%S")).days

    years = np.array([nc.num2date(t, units=time_var.units, calendar=getattr(time_var, 'calendar', 'gregorian')).year for t in time_var])

    unique_years = np.unique(years)

    mean_eto_per_year = np.zeros(len(unique_years))

    for i, year in enumerate(unique_years):
        mask = (years == year)
        mean_eto_per_year[i] = np.nanmean(eto_var[mask])
        
    if 'with_FAO-PM_eq' in file_path:
        co2_input = "without_CO2"
    elif 'with_CO2' in file_path:
       co2_input = "with_CO2"
    # Create a DataFrame for the current file
    file_df = pd.DataFrame({'GCM': gcm_name,'CO2':co2_input,
                            'Year': unique_years,
                            'Mean Ref. ETo': mean_eto_per_year})

    # Append the DataFrame to the result DataFrame
    result_df = pd.concat([result_df, file_df], ignore_index=True)

    dataset.close()

print("start_year:", start_year)
# Export the final DataFrame to Excel
excel_file_path = os.path.join(directory_path_var_2, start_year)
excel_file_path = excel_file_path+'_' + str(int(start_year) + 9) + 'annual_avg.xlsx'
# excel_file_path = 'E:\ISIMIP Climate Data\eto_masked_changed\difference_with_headings.xlsx'        # enter the file path where data has to be saved
result_df.to_excel(excel_file_path, index=False)

print(f"Data exported to {excel_file_path}")
