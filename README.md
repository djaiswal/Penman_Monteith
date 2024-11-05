# Reference Evapotranspiration (ETo) Calculation and Analysis

This project contains Python scripts for calculating and analyzing reference evapotranspiration (ETo) from NetCDF files. The scripts include codes for calculating ETo with and without considering the effect of CO₂ concentration changes, analyzing the relative and percentage changes over time, and generating plots for data visualization.

## Table of Contents

- [I. Setup Instructions](#i-setup-instructions)
- [II. Calculate Daily Net_radiation](#ii-calculate-daily-net-radiation)
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

## II. Calculate Daily Net radiation
`zenith_3652days.py`  
This script is used to replicate values of elevation for 3652 time steps (10 years). The new values are saved in new NetCDF file named elevation_data_for_3652_days.nc.

### Steps:
1. Select the netcdf file with the data for elevation for one time step (The code assumes that this file has the data of elevation in the variable name 'elevation').
2. Choose where to save the output file with the data of radiation for 3652 days.

   
`calculate_radiation.py`  
This script calculates daily radiation input NetCDF files. The calculated values are saved in new NetCDF files.

### Steps:
1. Select the required input files.
2. Choose where to save the output daily radiation files.
   
## III. Calculate Daily Reference Evapotranspiration

`calculate_Eto.py`  
This script calculates daily ETo using input NetCDF files. The calculated values are saved in new NetCDF files.

### Steps:
1. Create a directory to save the output files, e.g., `daily_ETo`.
2. Select the directory containing input climate data in NetCDF format.
3. Choose where to save the output daily ETo files.

## IV. Calculate Annual Values of Reference Evapotranspiration

### `annual_avg.py`
Calculates the annual average of daily ETo over a selected region and saves the results as an Excel file.

#### Steps:
1. Create a directory to save the Excel file, e.g., `annual_avg`.
2. Select the folder containing NetCDF files with daily ETo data.
3. Specify where to save the Excel file with the annual averages.

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

### `plot_4ETosubplots.py`
Generates subplots showing ETo with and without CO₂ effect for 2021-2030 and 2091-2100.

#### Steps:
1. Select NetCDF files for each ETo scenario for the chosen GCM.
2. Provide the shapefile of India and specify output directory for the plot.

### Plot Relative Change in Reference Evapotranspiration
Creates a plot showing the temporal mean of relative change in ETo for a single GCM over ten years.

#### Steps:
1. Select NetCDF files with relative change data for 2021-2030 and 2091-2100.
2. Provide the shapefile of India and specify the output directory.

### Plot Mean of Relative Change in Reference Evapotranspiration for Five GCMs
Plots the mean of relative change in ETo across all five GCMs.

#### Steps:
1. Select NetCDF files for each GCM and period.
2. Provide the shapefile of India and specify the output directory.
