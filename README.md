# Modified FAO-PM Equation

## Data collection

The datasets for the calculation of daily reference ETo were downloaded from
the ISIMIP repository for weather data. These datasets cover the CMIP6-based and bias-adjusted atmospheric climate input data of five CMIP6 GCMS, viz, GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0 and UKESM1-0-LL. These datasets belonging to each GCMs, if saved in separate folders, can serve better in the calculation and analysis for this study.

The data was downloaded as NetCDF files bound by latitudes 7.25 and 37.25 and longitudes 67.75 and 97.75, for the periods from 2021 to 2030 and from 2091 to 2100.

Note: The file name of the daily weather data is expected to include the name of the GCM and and its corresponding starting year.

The shapefile of India that was used to create plots and mask file for the NetCDF files downloaded as stated above was obtained from the site.


## Calculation of Reference Evapotranspiration

### calculate_ETo.py
	
This code calculates the reference evapotranspiration without considering the effect of CO2 and considering the effect of  CO2 as two NetCDF files respectively in a selected directory.

#### Steps 
Create a directory to save the NetCDF files with data of reference evapotranspiration. ( For Eg: Daily ETo NetCDF files)
Enter the starting year of the decade for which reference evapotranspiration has to be calculated. For eg. , if the reference evapotranspiration for 2021 to 2030 has to be calculated, enter 2021 in the terminal.
From the first file selection window, select the NetCDF files with the data for 
	i . Relative humidity
	ii.  Pressure
	iii.  Long radiation
	iv.  Short radiation
	v.  Surface windspeed
	vi. Daily Maximum Temperature
	vii. Daily Minimum Temperature
Click on the ‘Finish’ button to finish selecting the files.
From the second selection window, select the directory where the NetCDF files with reference evapotranspiration calculated without considering the effect of CO2 and calculated considering the effect of CO2 have to be saved ( folder named ‘Daily ETo NetCDF files’  as described in step II.a.1).

## Mask NetCDF files to project only the area inside the Indian Mainland.

A NetCDF file that can be used to mask the area inside the Indian Mainland,  can be created using the shapefile of India.

Using ArcGIS software, the shapefile of India can be processed to create a NetCDF file with value = 0 for the area inside India and value = 1 for the area outside. This is a NetCDF file with one timestep.

### mask0and1flipping.py 

This code creates a NetCDF file, that has value = 1 for the area inside India and values set as nan (not a number) for the area outside India using the NetCDF file obtained using ArcGIS software as mentioned above. The file is saved as ‘0and1proper_mask.nc’ in a selected folder. 

#### Steps:
Create a directory where the mask files have to be stored (E.g.: Mask files).
Select the NetCDF file obtained using ArcGIS software (mentioned earlier) from the first file selection window.
Select the directory where the NetCDF file which can be used to mask other NetCDF files is to be saved from the next selection window (the directory created as per step III.a.1). 

### mask3652days.py

The NetCDF mask file created using mask0and1flipping.py can mask only NetCDF files with only one timestep appropriately. To mask a NetCDF file with 3652 timesteps (10 years in this case), the values in the former NetCDF file are replicated for 3652 timesteps and saved as a NetCDF file.

#### Steps: 
Select the mask file for one timestep (mask file created using mask0and1flipping.py) from the first selection window.
Select the directory where the new mask file for 3652 timesteps (10 years) has to be saved (the directory created as per step III.a.1).

### mask_files_in_a_directory.py

The NetCDF files containing the reference evapotranspiration data can be masked appropriately to retain only the data of the area inside the Indian mainland using this code.

#### Steps:
Select the directory containing the NetCDF files of reference evapotranspiration to be masked and the mask NetCDF file from the selection window.


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



