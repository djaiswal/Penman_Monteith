
# NOTE : This code calculates annual avg of ETo from daily ETo for both the timeperiods and five GCMs.
#        The calculated values are stored in an excel sheet. 

#        Select the directory containing NetCDF files of daily Eto. from selection window.
#        The second input required is the absolute filepath to excel sheet to which 
#        annual average of Eto has to be saved. 


import netCDF4 as nc
import numpy as np
from datetime import datetime
import pandas as pd
import os
import glob
import tkinter as tk
from tkinter import filedialog

directory_path = ""
def browse_directory():
    global directory_path
    directory_path = filedialog.askdirectory()
    if directory_path:
        directory_label.config(text=f"Selected Directory with ETo files : {directory_path}")

def finish():
    root.destroy()

# Create the main window
root = tk.Tk()
root.title("Directory Browser to find annual average of relative change")
root.geometry("800x200")  # Set the window size to 800x200
directory_button = tk.Button(root, text="Browse Directory ", command=browse_directory)
directory_button.pack(pady=10)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=20)
finish_button = tk.Button(root, text="Finish selecting directory", command=finish)
finish_button.pack(pady=20)
root.mainloop()


# Define the pattern to match NetCDF files
file_pattern = '*.nc'

# Create a list of file paths
file_paths = glob.glob(os.path.join(directory_path, file_pattern))

# Initialize an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['File', 'Year', 'Mean Ref. ETo'])

# Iterate through each file
for file_path in file_paths:
    dataset = nc.Dataset(file_path, 'r')

    time_var = dataset.variables['time']
    eto_var = dataset.variables['evapotranspiration'][:]

    # Assuming time is in a standard datetime format, you can use the following to extract the year
    if '2021' in file_path:
        base_date = "2021-01-01 00:00:00"
    elif '2091' in file_path:
        base_date = "2091-01-01 00:00:00"
    base_datetime = nc.num2date(0, units=time_var.units, calendar=getattr(time_var, 'calendar', 'gregorian'))
    days_since_base_date = (base_datetime - datetime.strptime(base_date, "%Y-%m-%d %H:%M:%S")).days

    years = np.array([nc.num2date(t, units=time_var.units, calendar=getattr(time_var, 'calendar', 'gregorian')).year for t in time_var])

    unique_years = np.unique(years)

    mean_eto_per_year = np.zeros(len(unique_years))

    for i, year in enumerate(unique_years):
        mask = (years == year)
        mean_eto_per_year[i] = np.nanmean(eto_var[mask])

    # Create a DataFrame for the current file
    file_df = pd.DataFrame({'File': [os.path.basename(file_path)] * len(unique_years),
                            'Year': unique_years,
                            'Mean Ref. ETo': mean_eto_per_year})

    # Append the DataFrame to the result DataFrame
    result_df = pd.concat([result_df, file_df], ignore_index=True)

    dataset.close()

# Export the final DataFrame to Excel
excel_file_path = input("Enter the absolute filepath to the excel sheet where the data has to be saved : ")
# excel_file_path = 'E:\ISIMIP Climate Data\eto_masked_changed\difference_with_headings.xlsx'        # enter the file path where data has to be saved
result_df.to_excel(excel_file_path, index=False)

print(f"Data exported to {excel_file_path}")
