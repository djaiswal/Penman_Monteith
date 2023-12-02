import netCDF4 as nc

input_file = r"E:\ISIMIP Climate Data\SHapefiles and netcdf file_02_12_2023\0and1_proper_india_02_12_2023.nc"
data = nc.Dataset(input_file , 'r')

print("This file has" + str(len(data.variables.keys())) + "variables.")
print("The following are the variables in this file:")
print(data.variables.keys())

for i in range(len(data.variables.keys())):

    print("The following are the values of " + str(list(data.variables.keys())[i]))
    print(data.variables[str(list(data.variables.keys())[i])][:])

