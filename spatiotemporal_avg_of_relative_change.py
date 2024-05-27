# NOTE : This code calculates the average of difference in ETo over space and time.
#        The window to select the Netcdf files containing data of difference in time
#        appears on running the code. It is required to select the correct file for
#        each GCM to calculate the average accurately. After selecting all the files,
#        click on "Finish selecting files" button to close the window to complete the
#        calculation.

from netCDF4 import Dataset
import numpy as np
import os
import glob
import tkinter as tk
from tkinter import filedialog

filelist = []

# Assuming that we are calculating the spatio-temporal average
file_paths = [""] * 5
# of relative change for 5 GCMs together. If not, change the number here
# add or remove the GCM's name in the following line.
GCM_list = ['GFDL-ESM4', 'IPSL-CM6A-LR',
            'MRI-ESM2-0', 'MPI-ESM1-2-HR', 'UKESM1-0-LL']



# The following section of this module is to create an interface to select the files.

def browse_file(file_index):
    file_path = filedialog.askopenfilename()
    if file_path:
        file_paths[file_index] = file_path
        labels[file_index].config(
            text=f" {GCM_list[file_index]} : {file_path}")

def finish():
    root.destroy()

# Create the main window for selecting the NetCDF files.
root = tk.Tk()
root.title("File Browser")
root.geometry("800x400")  # Set the window size to 800x300

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



# The following part of the code is to calculate the spatiotemporal mean of ETo for each selected file.

print("\n\n")
print("Mean of difference in ETo with and without CO2 averaged over space and time")

j = 0

for file_path in file_paths:
    dataset = Dataset(file_path, 'r')

    eto_var = dataset.variables['difference_in_eto'][:]
    mean_difference = np.nanmean(eto_var, axis=(0, 1, 2))
    print(GCM_list[j], ' - ', mean_difference)
    if (j > 4):
        break
    j += 1
