from netCDF4 import Dataset
import numpy as np
import os
import glob

directory_path = 'E:\ISIMIP Climate Data\eto_masked_changed\difference_netcdf_files'
file_pattern = '*.nc'
file_paths = glob.glob(os.path.join(directory_path, file_pattern))

for file_path in file_paths:
    dataset = Dataset(file_path, 'r')

    eto_var = dataset.variables['difference_in_eto'][:]
    mean_difference = np.nanmin(eto_var, axis=(0,1,2))
    print("Mean of difference in ETo with and without CO2 averaged over space and time")
    print(file_path[73:file_path.find("_diff")], ' - ',{mean_difference} )
