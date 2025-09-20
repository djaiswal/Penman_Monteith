# NOTE : This code calculates the mean, maximum and minimum percentage change of ETo values
#        For this we need NetCDF files with data of percentage change of ETo for different GCMs under different decades.
from netCDF4 import Dataset
import numpy as np
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import os

# Assuming that we are calculating the spatio-temporal average of percentage change
# of percentage change for 5 GCMs together. If not, change the number here and
# add or remove the GCM's name in the following lines accordingly.
file_paths = [""] * 15
GCM_list = ['GFDL-ESM4 2021-2030','GFDL-ESM4 2051-2060','GFDL-ESM4 2091-2100',
            'IPSL-CM6A-LR 2021-2030','IPSL-CM6A-LR 2051-2060','IPSL-CM6A-LR 2091-2100',
            'MPI-ESM1-2-HR 2021-2030 ', 'MPI-ESM1-2-HR 2051-2060', 'MPI-ESM1-2-HR 2091-2100', 
            'MRI-ESM2-0 2021-2030 ','MRI-ESM2-0 2051-2060','MRI-ESM2-0 2091-2100 ',
            'UKESM1-0-LL 2021-2030 ', 'UKESM1-0-LL 2051-2060', 'UKESM1-0-LL 2091-2100 ']

# The following section of this module is to create an interface to select the files.

def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {GCM_list[file_index]} : {file_path}")
def browse_directory():
    global directory_path_var_1
    path = filedialog.askdirectory()
    if path:
        directory_path_var_1.set(path)
        directory_label.config(text=f"Selected Directory with ETo files : {directory_path_var_1.get()}")

def finish():
    root.destroy()

# Create the main window for selecting the NetCDF files.

root = tk.Tk()
root.title("File Browser")
root.geometry("800x750")  # Set the window size to 800x300
directory_path_var_1 = tk.StringVar()
labels = []
buttons = []

for i in range(len(file_paths)):

    button = tk.Button(
    root, text=f"Browse File for {GCM_list[i]}", command=lambda i=i: browse_file(i))
    button.pack(pady=5)
    buttons.append(button)

    label = tk.Label(root, text=f"{GCM_list[i]}: ")
    label.pack(pady=5)
    labels.append(label)
    

    # Create a button to finish selecting files and close the window
finish_button = tk.Button(root, text="Finish Selecting Files", command=finish)
finish_button.pack(pady=5)
root.mainloop()

root = tk.Tk()
root.title("File Browser")
root.geometry("800x200")
directory_button = tk.Button(root, text="Browse Directory to save spatiotemporal mean of percentage change in ETo", command=browse_directory)
directory_button.pack(pady=5)
directory_label = tk.Label(root, text="Selected Directory: ")
directory_label.pack(pady=5)

finish_button = tk.Button(root, text="Finish Selecting Files", command=finish)
finish_button.pack(pady=5)
root.mainloop()
directory_path_var_1 = directory_path_var_1.get()


# The following part of the code is to calculate the spatiotemporal mean of ETo for each selected file.

result_df = pd.DataFrame(columns=['GCM', 'Years', 'Spatiotemporal Mean Ref. ETo percentage change', 'Max during the decade', 'Min during the decade', 'Max in the temporal mean', 'Min in temporal mean'])

mean_difference = []

for file_path in file_paths:
    dataset = Dataset(file_path, 'r')
    if '2021' in file_path:
        start_year = "2021-2030"
    elif '2051' in file_path:
        start_year = "2051-2060"
    elif '2091' in file_path:
        start_year = "2091-2100"

    if ('gfdl' in file_path) or 'GFDL' in file_path:
        gcm_name = 'GFDL-ESM4'
    elif ('ipsl' in file_path) or 'IPSL' in file_path:
        gcm_name = 'IPSL-CM6A-LR'
    elif 'mpi' in file_path or 'MPI'in file_path:
        gcm_name = 'MPI-ESM1-2-HR'
    elif 'mri' in file_path or 'MRI'in file_path:
        gcm_name = 'MRI-ESM2-0'
    elif 'ukesm' in file_path or 'UKESM'in file_path:
        gcm_name = 'UKESM1-0-LL'


    eto_var = dataset.variables['percentge_change'][:]
    file_df = pd.DataFrame({'GCM': gcm_name,
                            'Years': start_year,
                            'Spatiotemporal Mean Ref. ETo percentage change': [np.nanmean(eto_var, axis=(0, 1, 2))],
                            'Max during the decade' : [np.nanmax(eto_var)],
                            'Min during the decade' : [np.nanmin(eto_var)], 
                            'Max in the temporal mean' : [np.nanmax(np.nanmean(eto_var, axis=(0)))],
                            'Min in temporal mean':[np.nanmin(np.nanmean(eto_var, axis=(0)))]})
    
    # Append the DataFrame to the result DataFrame
    result_df = pd.concat([result_df, file_df], ignore_index=True)

    dataset.close()

excel_file_path = os.path.join(directory_path_var_1, 'spatiotemporal_avg_percentage_change.xlsx')
result_df.to_excel(excel_file_path, index=False)

print(f"Data exported to {excel_file_path}")