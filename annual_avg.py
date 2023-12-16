import os
import glob
import netCDF4 as nc
import numpy as np
from datetime import datetime
import pandas as pd

# Directory containing NetCDF files
directory_path = 'E:\ISIMIP Climate Data\eto_masked_changed\eto_files'

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
excel_file_path = 'E:\ISIMIP Climate Data\eto_masked_changed\difference_with_headings.xlsx'
result_df.to_excel(excel_file_path, index=False)

print(f"Data exported to {excel_file_path}")
