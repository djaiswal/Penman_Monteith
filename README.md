# Reference Evapotranspiration (ETo) Calculation and Analysis

This project contains Python scripts for calculating and analyzing reference evapotranspiration (ETo) from NetCDF files. The scripts include codes for calculating ETo with and without considering the effect of CO₂ concentration changes, analyzing the relative and percentage changes over time, and generating plots for data visualization.

## Table of Contents

- [I. Setup Instructions](#i-setup-instructions)
- [II. Calculate Daily Net_radiation](#ii-calculate-daily-net-radiation) (Optional : Do this if only net radiation files have to be generated)
- [III. Calculate Daily Reference Evapotranspiration](#iii-calculate-daily-reference-evapotranspiration)
- [IV. Calculate Annual Values of Reference Evapotranspiration](#iv-calculate-annual-values-of-reference-evapotranspiration)
- [V. Calculate the Relative Change in Reference Evapotranspiration](#v-calculate-the-relative-change-in-reference-evapotranspiration)
- [VI. Calculate Reference Evapotranspiration for the Four Seasons](#vi-calculate-reference-evapotranspiration-for-the-four-seasons)
- [VII. Calculate Percentage Change in Reference Evapotranspiration](#vii-calculate-percentage-change-in-reference-evapotranspiration)
- [VIII. Plotting the Data](#viii-plotting-the-data)

---

## I. Setup Instructions

1. Install necessary Python packages:
    ```bash
    pip install netCDF4 xarray pandas numpy matplotlib
    ```
2. Download NetCDF files and set up directories as described in each section.

## II. Calculate Daily Net radiation (Optional : Do this if only net radiation files have to be generated)
`calculate_elevation.py`  
This script is used to replicate values of elevation for time steps specified by the user. The new values are saved in new NetCDF file named elevation_data_for_<timesteps>_days.nc. 


### Steps:
1. Run the code `calculate_elevation.py`. Inorder to run this code successfully, you need to have a file with the data of elevation for atleast one time step. This file will provide the input elevation data. The code assumes that this file has the data of elevation in the variable name `'elevation'`.
2.  An interface appears to select the netcdf file with input elevation data and the number of timesteps for replicating the data. 
3. Select the netcdf file with the data for elevation for one time step.
4. Enter the number of time steps for which the data has to be replicated. In this project, we calculate the data for ten consecutive years in most cases. So here the expected value is the number of days in ten years.i.e, 3652 for any consecutive ten years not including the year 2100, or 3653 for consecutive ten years including 2100. Click on `Finish selecting file` button.
   <img width="1001" height="348" alt="image" src="https://github.com/user-attachments/assets/15708dbf-9cd8-4bfc-9a5d-00b728716e3d" />

5. From the next interface, choose where to save the output file with the data of radiation for the given number of days. The netcdf file with the data of elevation for the specified number of days will be saved in this folder with the name elevation_data_for_<timesteps>_days.nc.

   <img width="997" height="282" alt="image" src="https://github.com/user-attachments/assets/6f4182fa-876e-4bc0-b46f-e2c64ecc9718" />


   
`calculate_radiation.py`  
This script calculates daily radiation input NetCDF files. The calculated values are saved in new NetCDF files.

### Steps:
1. Select the required input files.
2. Choose where to save the output daily radiation files.
   
## III. Calculate Daily Reference Evapotranspiration

`calculate_Eto_ssp126.py`  
This script calculates daily ETo using input NetCDF files for the scenario SSP-1-2.6. The calculated values are saved in new NetCDF files.

### Steps:
1. Create a directory to save the output files, e.g., `ETo_files`.
2. Run the code. 
3. Enter the starting year of the decade for which reference evapotranspiration has to be calculated.
   i.e, 2021 for the decade 2021-2030, 2051 for the decade 2051-2060, or 2091 for the decade 2091-2100.
4. Select the files containing the weather data (for SSP1-2.6) from the first window. 
5. Select the parent folder where the ETo files for each GCM has to be saved. (For Eg: ETo_files)
        Separate folders are created for different GCMs and respective files are saved there.
        This creates a directory structure as follows:
 
           ETo_files
                |--- GFDL-ESM4 (contains ETo files for GFDL-ESM4)
                |--- IPSL-CM6A-LR (contains ETo files for IPSL-CM6A-LR)
                |--- MPI-ESM1-2-HR (contains ETo files for MPI-ESM1-2-HR)
                |--- MRI-ESM2-0 (contains ETo files for MRI-ESM2-0)
                |--- UKESM1-0-LL (contains ETo files for UKESM1-0-LL)

## IV. Calculate Annual Values of Reference Evapotranspiration

### `annual_avg.py`
Calculates the annual average of daily ETo over a selected region and saves the results as an Excel file.

#### Steps:
1. Select the directory containing NetCDF files of daily Eto from selection window.
        Suppose the directory structure is as follows:
           ETo_files
                |--- GFDL-ESM4
                |--- IPSL-CM6A-LR
                |--- MPI-ESM1-2-HR
                |--- MRI-ESM2-04-0
                |--- UKESM1-0-LL

        Select the directory 'ETo_files' in the selection window.

2.        The second input required is the directory to which annual average of Eto has to be saved.

### `annual_sum_of_eto.py`
Calculates the total annual ETo for each NetCDF file, storing results in an Excel sheet.

#### Steps:
1. Create a directory to save the output, e.g., `annual_sum`.
2. Select the folder with daily ETo NetCDF files.
3. Specify where to save the Excel files.

## V. Calculate the Relative Change in Reference Evapotranspiration

### `relative_change_eto_gcm.py`
Calculates the relative change in ETo due to CO₂ concentration changes, saving output as NetCDF files for 2021-2030 and 2091-2100.

#### Steps:
1. Create a directory to save the output NetCDF files, e.g., `daily_relative_change_in_ETo`.
2. Select the ETo NetCDF files:
   - With and without CO₂ effect for 2021-2030
   - With and without CO₂ effect for 2091-2100

### `spatiotemporal_avg_of_relative_change.py`
Calculates the spatiotemporal average of the relative change in ETo across all 5 GCMs for 2021-2030 and 2091-2100. Results are printed on the terminal.

#### Steps:
1. Create a directory to save the Excel output, e.g., `spatiotemporal_avg`.
2. Select the NetCDF files for each GCM for the specified years.

### `annual_avg_of_relative_change.py`
Calculates the annual average of the relative change in ETo and saves it as an Excel sheet.

#### Steps:
1. Select the folder containing NetCDF files with relative change data.

## VI. Calculate Reference Evapotranspiration for the Four Seasons

`seasonal_avg.py`  
This code calculates seasonal ETo averages for each GCM for specified years, saving data as NetCDF and Excel files, and generating plots.

### Steps:
1. Create a directory, e.g., `seasonal_ETo`, with subdirectories:
   - `seasonal_netcdf_files`
   - `seasonal_plots`
   - `seasonal_spatiotemporal_avg`
2. Select the appropriate NetCDF files (with and without CO₂ effect) for:
   - 2021-2030
   - 2091-2100
3. Provide the shapefile of India in .shp format for visualization.

## VII. Calculate Percentage Change in Reference Evapotranspiration

### `percentage_change_netcdf.py`
Calculates the percentage change in ETo due to CO₂ effect.

#### Steps:
1. Select NetCDF files with ETo data for 2021-2030 and 2091-2100 (both with and without CO₂ effect).
2. Specify where to save the NetCDF output with percentage change data.

### `percentage_change_min_max_mean.py`
Calculates:
- Spatiotemporal mean of percentage change
- Maximum and minimum values over ten years

#### Steps:
1. Select files with percentage change data for specified years across GCMs.
2. Specify where to save the Excel output.

## VIII. Plotting the Data

Plotting scripts are organized in the `code_for_plots` directory by figure and scenario. Use the following scripts for generating publication-quality plots:

### Fig4: Decadal Average of ETo
Script:
- `both_scenario_gcm_v1.py` (comparison across scenarios)
Location: `code_for_plots/Fig4-Decadal_avg_of_ETo/`

### Fig5: Decadal Average of Difference in ETo
Script:
- `plot_both_gcms_v2.py` (plots both GCMs)
Location: `code_for_plots/Fig5-Decadal_avg_of_difference_in_ETo/`

### Fig6: Annual Time Series of ETo
Scripts:
- `plot_both_gcms.py` (plots both GCMs)
- `plotssp126.py` (SSP1-2.6 scenario)
- `plotssp585.py` (SSP5-8.5 scenario)
Location: `code_for_plots/Fig6-Annual_time_series_of_ETo/`

### Fig7: Seasonal Average of Relative Change in ETo
Script:
- `plot_both_gcms_v1.py` (seasonal relative change)
Location: `code_for_plots/Fig7-Seasonal_average_of_relative_change_in_ETo/`

### Fig9: Spatiotemporal Average of ETo
Currently, this folder is empty. Add scripts here for spatiotemporal plots as needed.

#### General Steps for Plotting:
1. Select the appropriate NetCDF files for the scenario, GCM, and period as required by the script.
2. Provide the shapefile of India (see `code_for_plots/Shapefiles/`).
3. Specify the output directory for saving plots.
4. Run the desired plotting script from its respective folder.

Refer to comments in each script for specific input requirements and customization options.
