# Reference Evapotranspiration Calculation

This repository contains scripts for calculating daily reference evapotranspiration (ETo) from climate data, specifically tailored for data within the Indian mainland. The project utilizes datasets from the ISIMIP repository and processes NetCDF files to produce masked data relevant to the region.

## Table of Contents
1. [Data Collection](#data-collection)
2. [Calculation of Reference Evapotranspiration](#calculation-of-reference-evapotranspiration)
3. [Masking NetCDF Files](#masking-netcdf-files-to-project-only-the-area-inside-the-indian-mainland)

## Data Collection

The datasets for calculating daily reference ETo were downloaded from the ISIMIP repository, covering bias-adjusted atmospheric climate data from five CMIP6 GCMs:
- GFDL-ESM4
- IPSL-CM6A-LR
- MPI-ESM1-2-HR
- MRI-ESM2-0
- UKESM1-0-LL

Each dataset is stored in separate folders for easier analysis. Data files are in NetCDF format, spanning latitudes 7.25 to 37.25 and longitudes 67.75 to 97.75, covering two periods: 2021-2030 and 2091-2100.

> **Note**: The filenames for daily weather data must include the GCM name and the starting year for proper processing.

The shapefile of India, used to create plots and mask files, was obtained from an external source.

## Calculation of Reference Evapotranspiration

The `calculate_ETo.py` script calculates daily reference evapotranspiration with and without considering CO₂ effects. It outputs two NetCDF files for each dataset in the selected directory.

### Steps
1. **Create Output Directory**  
   Create a directory for saving the NetCDF files with the reference evapotranspiration data. (e.g., `Daily ETo NetCDF files`)
   
2. **Run the Script**  
   Run the script and enter the starting year of the target decade (e.g., 2021 for 2021-2030).
   
3. **Select Input Files**  
   Select the following NetCDF files in the prompted file selection windows:
   - Relative Humidity
   - Pressure
   - Longwave Radiation
   - Shortwave Radiation
   - Surface Windspeed
   - Daily Maximum Temperature
   - Daily Minimum Temperature

4. **Select Output Directory**  
   Select the directory to save the calculated ETo NetCDF files for CO₂ and non-CO₂ calculations.

## Masking NetCDF Files to Project Only the Area Inside the Indian Mainland

A mask for the Indian mainland is created using a shapefile, with steps provided to ensure accurate masking of multi-timestep NetCDF files.

### Creating the Mask File

Using ArcGIS software, process the Indian shapefile to produce a NetCDF file with a value of `0` inside India and `1` outside.

### Scripts and Usage

1. **`mask0and1flipping.py`**  
   Converts the ArcGIS mask file to a format suitable for masking other NetCDF files, setting values inside India to `1` and outside to `NaN`.

   **Steps**:
   - Create a directory to save mask files (e.g., `Mask files`).
   - Select the original mask NetCDF file created with ArcGIS.
   - Select the output directory for the processed mask file (`0and1proper_mask.nc`).

2. **`mask3652days.py`**  
   Converts a single-timestep mask to cover 3652 timesteps, suitable for a 10-year period.

   **Steps**:
   - Select the one-timestep mask file created by `mask0and1flipping.py`.
   - Select the output directory for the 3652-timestep mask file.

3. **`mask_files_in_a_directory.py`**  
   This script masks all NetCDF files containing ETo data, retaining only data within the Indian mainland.

   **Steps**:
   - Select the directory containing ETo NetCDF files.
   - Select the mask NetCDF file to apply to each file in the directory.


## IV. Calculate annual values of reference evapotranspiration.

## annual_avg.py

This code calculates the annual average of reference evapotranspiration over space for all NetCDF files with daily data of reference evapotranspiration in a selected directory and saves it as an Excel file.
	
####	Steps :
	
Create a directory to save the output Excel files (For Eg: annual avg).
Select the directory containing the NetCDF files with daily data of reference evapotranspiration. 
Select the directory where the Excel file with the annual average of reference evapotranspiration data has to be saved (directory created as per step IV.a.1).

### annual_sum_of_eto.py

This code calculates the total reference evapotranspiration for every year for each NetCDF file (with the daily data of reference evapotranspiration)  in a selected folder. The calculated values are stored in an Excel sheet.

#### Steps :
Create a folder to save the Excel files with the data of the annual sum of reference evapotranspiration (For Eg: annual sum)
Select the folder that contains the NetCDF files with the daily data of reference evapotranspiration 
Select the folder where the Excel files are to be saved from the directory selection window (folder created as per step IV.b.1). 


## V. Calculate the relative change in reference evapotranspiration.
	
Optional step : Create a directory to save all the files related to relative change in reference evapotranspiration.Create three subdirectories in this folder as mentioned in steps V.a.1 and V.b.1

### relative_change_eto_gcm.py

This code calculates the relative change between the evapotranspiration calculated without considering the effect of CO2 and considering the effect of CO2 for a GCM. The datasets (for 2021-2030 and 2091-2100) that have this relative change are saved as two NetCDF files respectively in a selected folder.

#### Steps : 
Create a directory where the NetCDF files with the relative change in reference evapotranspiration are to be saved. (For Eg : daily relative change in ETo). 
Select the ETo NetCDF files with the value of reference evapotranspiration calculated 
	i.  considering the effect of CO2 for 2021-2030
	ii. without considering the effect of CO2 for 2021-2030
	iii. considering the effect of CO2 for 2091-2100
iv. without considering the effect of CO2 for 2091-2100
v. the directory where the NetCDF files with the relative change in reference evapotranspiration are to be saved.

### spatiotemporal_avg_of_relative_change.py
	
This code calculates the average of relative change in ETo over space and time for all 5 GCMs for the years 2021-2030 and 2091-2100. These values are printed on the terminal.


#### Steps : 
Create a subdirectory where the Excel file with the data of spatiotemporal average has to be saved. 
From the file selection window, select the NetCDF files with data of relative change in reference evapotranspiration for 
i. GFDL-ESM4 2021-2030
ii. GFDL-ESM4 2091-2100
iii. IPSL-CM6A-LR 2021-2030
iv. IPSL-CM6A-LR 2091-2100
v. MPI-ESM1-2-HR 2021-2030
vi. MPI-ESM1-2-HR 2091-2100
vii. MRI-ESM2-0 2021-2030
viii. MRI-ESM2-0 2091-2100
ix. UKESM1-0-LL 2021-2030
x. UKESM1-0-LL 2091-2100
xi. directory where the Excel file with the data of spatiotemporal average has to be saved (created in step V.b.1)

### annual_avg_of_relative_change.py

This code calculates the annual average of relative change in reference evapotranspiration (ETo calculated with and without the effect of CO2) and saves the data as an Excel sheet.

#### Steps : 
Select the folder that contains the NetCDF files with data of relative change in reference evapotranspiration (these files can be created by running relative_change_eto_gcm.py) 

## V. Calculate reference evapotranspiration for the four seasons
### seasonal_avg.py
This code calculates the reference evapotranspiration for the four seasons from the NetCDF files containing the daily data of reference evapotranspiration for ten years for a GCM.
The daily data of evapotranspiration for four different seasons are extracted and are saved as four different NetCDF files, for each NetCDF file selected (NetCDF file with data of daily reference evapotranspiration for ten years) in a subdirectory named seasonal_netcdf_files.
 The spatiotemporal mean of reference evapotranspiration for each season for the selected NetCDF files are calculated and saved in an Excel file in a subdirectory.
 The plots showing the temporal mean of reference evapotranspiration for four different seasons are also saved in a different subdirectory 


Create a directory for saving the files of reference evapotranspiration for different seasons, say ‘seasonal_ETo’.  In this parent directory, three subdirectories,
	i. seasonal_netcdf_files
	ii. seasonal_plots 
	iii. seasonal_spatiotemporal_avg
will be created after running the code to save the corresponding files.

#### Steps :
After running seasonal_avg.py, select
i. The  NetCDF file with the data of ETo calculated without considering the effect      of  CO2 for 2021-2030 
ii. The NetCDF file with the data of ETo calculated considering the effect of CO2 for 2021-2030
iii. The  NetCDF file with the data of ETo calculated without considering the effect of  CO2 for 2091-2100
iv. The NetCDF file with the data of ETo calculated considering the effect of CO2 for 2091-2100
v. the parent directory where the output files have to be saved
vi. The shapefile of India in .shp format

Note: The four selected NetCDF files should correspond to the same GCM.

## VI. Calculate the percentage of change in reference evapotranspiration

### percentage_change_netcdf.py

This code calculates the percentage of change in reference evapotranspiration calculated without considering the effect of CO2 and considering the effect of CO2.

#### Steps: 
Select the NetCDF files of  a GCM with the following data of reference evapotranspiration for calculated 
i. without considering the effect of CO2  for 2021-2030
ii. considering the effect of CO2 for 2021-2030 and 2091-2100.
iii.  without considering the effect of CO2  for 2091-2100
iv. considering the effect of CO2 for 2091-2100.

Then select the directory where the NetCDF file with the percentage of relative change has to be saved.
### percentage_change_min_max_mean.py

For every NetCDF file selected, this code calculates
i. The spatiotemporal mean of percentage change in reference evapotranspiration
ii. maximum percentage change of reference evapotranspiration during the entire ten years of consideration
iii.  minimum percentage change of reference evapotranspiration during the entire ten years of consideration
iv.  maximum percentage change of reference evapotranspiration from the temporal mean of percentage change in reference evapotranspiration (the mean over time of percentage change in ETo for ten years is taken first and the maximum value from this data is then found out)
v.  minimum percentage change of reference evapotranspiration from the temporal mean of percentage change in reference evapotranspiration (the mean over time of percentage change in ETo for ten years is taken first and the minimum value from this data is then found out)

####	Steps : 
Select the files with the percentage change in reference evapotranspiration for the years 2021-2030 and 2091-2100(for 5 GCMs) from the first selection window.
Select the directory where the Excel file with the output data has to be saved, from the second selection window. 

## VII. Plotting the data

### plot_4ETosubplots.py

This code creates a plot containing the following subplots with the temporal average of reference evapotranspiration calculated 
	i.  Without considering the effect of CO2 for 2021-2030
	ii. Considering the effect of CO2 for 2021-2030
	iii. Without considering the effect of CO2 for 2091-2100
	iv. Considering the effect of CO2 for 2091-2100



#### Steps : 
Select the following files belonging to a GCM with data of evapotranspiration calculated as follows : 
		i. Without considering the effect of CO2 for 2021-2030
		ii. Considering the effect of CO2 for 2021-2030
		iii. Without considering the effect of CO2 for 2091-2100
		iv. Considering the effect of CO2 for 2091-2100
		
Select the shapefile of India and then the directory where the plots have to be stored.

### Plot relative change in reference evapotranspiration
	 
This code plots the temporal mean of relative change in reference evapotranspiration (averaged over ten years) for a GCM

Select the NetCDF files with the relative change in reference evapotranspiration for 2021-2030 and 2091-2100. 
Select the shapefile of India and the directory where the plots have to be stored. 

Plot 	
This code plots the mean of relative change in reference evapotranspiration for five GCMs (relative change in reference ETo averaged over five GCMs)

Select the following NetCDF files with the data of relative change in reference evapotranspiration 
	GCM 			Year
	GFDL-ESM4  		2021-2030
GFDL-ESM4 		2091-2100
IPSL-CM6A-LR 	2021-2030
IPSL-CM6A-LR 	2091-2100
MPI-ESM1-2-HR 	2021-2030
MPI-ESM1-2-HR 	2091-2100
MRI-ESM2-0		2021-2030
MRI-ESM2-0 		2091-2100
UKESM1-0-LL 		2021-2030
UKESM1-0-LL 		2091-2100

Select the shapefile and the directory where the plot has to be saved. 



