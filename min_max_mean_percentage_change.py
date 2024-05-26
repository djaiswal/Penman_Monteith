from netCDF4 import Dataset
import numpy as np

withoutco2_2021 = Dataset(r"E:\ISIMIP Climate Data\et0_masked\masked_gfdl_et0_2021-2030_without_CO2.nc",'r')
withco2_2021 = Dataset(r"E:\ISIMIP Climate Data\et0_masked\masked_ipsl_et0_2091-2100_without_CO2.nc", 'r')
withoutco2_2091 = Dataset(r"E:\ISIMIP Climate Data\et0_masked\masked_ipsl_et0_2091-2100_without_CO2.nc",'r')
withco2_2091 = Dataset(r"E:\ISIMIP Climate Data\et0_masked\masked_ipsl_et0_2091-2100_with_CO2.nc", 'r')


withoutco2_data_2021 = withoutco2_2021.variables['evapotranspiration'][:]
withco2_data_2021 = withco2_2021.variables['evapotranspiration'][:]
withoutco2_data_2091 = withoutco2_2091.variables['evapotranspiration'][:]
withco2_data_2091 = withco2_2091.variables['evapotranspiration'][:] 

percentage_change_data_2021 = (withoutco2_data_2021-withco2_data_2021)/(withoutco2_data_2021)*100
percentage_change_data_2091 = (withoutco2_data_2091-withco2_data_2091)/(withoutco2_data_2091)*100

mean_percentchange_2021 = np.mean(percentage_change_data_2021, axis=(0,1,2))
mean_percentchange_2091 = np.mean(percentage_change_data_2091, axis=(0,1,2))

print("GCM : GFDL")
print('Percentage change in 2021-2030 is ', (mean_percentchange_2021))
print('Percentage change in 2091-2100 is ', (mean_percentchange_2091))
print('Max value (in 2021-2030) without CO2 is' , np.nanmax(np.mean(withoutco2_data_2021, axis=(0))))
print('Max value (in 2021-2030) with CO2 is' , np.nanmax(np.mean(withco2_data_2021, axis=(0))))
print('Max value (in 2091-2100) without CO2 is' , np.nanmax(np.mean(withoutco2_data_2091, axis=(0))))
print('Max value (in 2091-2100) with CO2 is' , np.nanmax(np.mean(withco2_data_2091, axis=(0))))


