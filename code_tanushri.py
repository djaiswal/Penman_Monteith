import math
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def fao56_penman_monteith_modified(net_rad, tmax, tmin, ws, rh_mean, atmos_pres, con_co2=422, shf=0.0, gs_ref=0.0061, co2_ref=330):
    """
    Estimate reference evapotranspiration (ETo) from a hypothetical
    short grass reference surface using the FAO-56 Penman-Monteith equation.
    Based on equation 6 in Allen et al (1998).
    :param net_rad: Net radiation at crop surface [MJ m-2 day-1]. If
        necessary this can be estimated using ``net_rad()``.
    :param t: Air temperature at 2 m height [deg Kelvin].
    :param ws: Wind speed at 2 m height [m s-1]. If not measured at 2m,
        convert using ``wind_speed_at_2m()``.
    :param svp: Saturation vapour pressure [kPa]. Can be estimated using
        ``svp_from_t()''.
    :param avp: Actual vapour pressure [kPa]. Can be estimated using a range
        of functions with names beginning with 'avp_from'.
    :param delta_svp: Slope of saturation vapour pressure curve [kPa degC-1].
        Can be estimated using ``delta_svp()``.
    :param psy: Psychrometric constant [kPa deg C]. Can be estimatred using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.

    :param con_co2: constration of co2 [ppm]. Can be estimatred using
        ``psy_const_of_psychrometer()`` or ``psy_const()``.

    :param gs : Conductivity of stoma, surface or stomatal conductance (m s-1)

    :param shf: Soil heat flux (G) [MJ m-2 day-1] (default is 0.0, which is
        reasonable for a daily or 10-day time steps). For monthly time steps
        *shf* can be estimated using ``monthly_soil_heat_flux()`` or
        ``monthly_soil_heat_flux2()``.
    :return: Reference evapotranspiration (ETo) from a hypothetical
        grass reference surface [mm day-1].
    :rtype: float
    
    """
    # calculating GS
    gs = gs_ref * (1 / (1 + 0.663 * ((con_co2 / co2_ref) - 1)))

    # Estimate mean daily temperature
    ''' Estimate mean daily temperature from the daily minimum and maximum temperatures.
    :param tmin: Minimum daily temperature [deg C]
    :param tmax: Maximum daily temperature [deg C]
    :return: Mean daily temperature [deg C]
    :rtype: float'''
    temperatures = (tmax + tmin) / 2.0 

    # Calculating Saturated Vapor Pressure (ES)
    svp = 0.6108 * np.exp((17.27 * temperatures) / (temperatures + 237.3))

    # Calculating Saturated Vapor Pressure at Tmin
    svp_tmin = 0.6108 * np.exp((17.27 * tmin) / (tmin + 237.3))

    # Calculating Saturated Vapor Pressure at Tmax
    svp_tmax = 0.6108 * np.exp((17.27 * tmax) / (tmax + 237.3))

    # Estimate actual vapour pressure (AVP)
    ''' Estimate actual vapour pressure (*ea*) from saturation vapour pressure at
    daily minimum and maximum temperature, and mean relative humidity.
    Based on FAO equation 19 in Allen et al (1998).
    :param svp_tmin: Saturation vapour pressure at daily minimum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param svp_tmax: Saturation vapour pressure at daily maximum temperature
        [kPa]. Can be estimated using ``svp_from_t()``.
    :param rh_mean: Mean relative humidity [%] (average of RH min and RH max).
    :return: Actual vapour pressure [kPa]
    :rtype: float'''
    avp = (rh_mean / 100.0) * ((svp_tmax + svp_tmin) / 2.0)

    # Calculate the psychrometric constant
    ''' Calculate the psychrometric constant.
    This method assumes that the air is saturated with water vapour at the
    minimum daily temperature. This assumption may not hold in arid areas.
    Based on equation 8, page 95 in Allen et al (1998).
    :param atmos_pres: Atmospheric pressure [kPa]. Can be estimated using
        ``atm_pressure()``.
    :return: Psychrometric constant [kPa degC-1].
    :rtype: float'''
    psy = 0.000665 * atmos_pres

    # Calculate es for each temperature
    es_values = []
    for temperature_min_val, temperature_max_val in zip(tmin, tmax):
        temperature = (temperature_min_val + temperature_max_val) / 2.0
        es = 0.6108 * np.exp((17.27 * temperature) / (temperature + 237.3))
        #es = round(es, 2)  # Round to two decimal points
        print(f"{temperature}={es}")
        es_values.append(es)

    # Create a line graph for es and temperature
    plt.plot(temperatures, es_values, label='Saturated Vapor Pressure')
    plt.xlabel('Temperature (Degree Celsius)')
    plt.ylabel('Saturated Vapor Pressure (kPa)')
    plt.title('Saturated Vapor Pressure vs Temperature')
    plt.legend()
    plt.grid(True)
    plt.show()
    

    # Calculate the slope using linear regression
    delta_svp,_ = np.polyfit(temperatures, es_values, 1)
    print("Slope:", delta_svp)

    # Modified Penman-Monteith equation
    a1 = (0.408 * (net_rad - shf) * delta_svp) / (delta_svp + (psy * ((1 + 0.035 * ws) / gs)))
    a2 = (900 * ws / temperatures) * (svp - avp) * psy / (delta_svp + (psy * ((1 + 0.035 * ws) / gs)))
    eto = a1 + a2

    return eto

# Read the netcdf file(1)
data = Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_hurs_ind_daily_2021_2030.nc')
data_1 = Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_ps_ind_daily_2021_2030.nc')
data_2 = Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rlds_ind_daily_2021_2030.nc')
data_3 =Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_rsds_ind_daily_2021_2030.nc')
data_4 = Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_sfcwind_ind_daily_2021_2030.nc')
data_5 = Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasmax_ind_daily_2021_2030.nc')
data_6 = Dataset(r'D:\IITPKD\Penman_Monteith\weatherData\gfdl-esm4_r1i1p1f1_w5e5_ssp585_tasmin_ind_daily_2021_2030.nc')



net_red_1 = data_2.variables['rlds'][:]
net_red_2 = data_3.variables['rsds'][:]

net_rad = net_red_1
'''# Reading The Excel Sheet Rh_Eto
Rh_Eto_data = pd.read_excel('C:\Python\Project_1\WeatherForPyeto.xlsx' , sheet_name='Rh_Eto')

# Reading The Excel Sheet Tmin_Eto
Tmin_Eto_data = pd.read_excel('C:\Python\Project_1\WeatherForPyeto.xlsx' , sheet_name='Tmin_Eto')

# Reading The Excel Sheet Tmax_Eto
Tmax_Eto_data = pd.read_excel('C:\Python\Project_1\WeatherForPyeto.xlsx' , sheet_name='Tmax_Eto')'''

print(data)
data.variables.keys()
lats = data.variables['lat'][1:10]
lons = data.variables['lon'][1:10]
time = data.variables['time'][1:10]

rh_mean = data.variables['hurs'][1:10]
atmos_pres = data_1.variables['ps'][1:10]


#net_rad = data_3.variables['air'][1:10]
ws = data_4.variables['sfcwind'][1:10]
temperature_min = data_6.variables['tasmin'][1:10]
temperature_max = data_5.variables['tasmax'][1:10]

#temperature_max = temperature_max-273.15
#temperature_min = temperature_min - 273.15

#con_co2 = data_3.variables['climatology_bounds'][:]


# Example usage of fao56_penman_monteith_modified function
eto = fao56_penman_monteith_modified(net_rad, temperature_max, temperature_min, ws, rh_mean, atmos_pres, )
print("Reference Evapotranspiration (ETo):", eto)

# Create a line graph for Co2 Concentration and Evapotranspiration 
'''plt.plot(con_co2, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Co2 Concentration (ppm)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Co2 Concentration vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()'''

# Create a line graph for Relative Humidity and Evapotranspiration 
plt.plot(rh_mean, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Relative Humidity (%)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Relative Humidity vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()

# Create a line graph for Tmin and Evapotranspiration 
plt.plot(temperature_min, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Minimum Temperature (Degree Celsius)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Minimum Temperature  vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()

# Create a line graph for Tmax and Evapotranspiration 
plt.plot(temperature_max, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Maximum Temperature (Degree Celsius)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Maximum Temperature  vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()


mp = Basemap(projection= 'merc',
             llcrnrlon= -50.0,
             llcrnrlat= 10.0,
             urcrnrlon= -95.0,
             urcrnrlat= 45.0,
             resolution='i')

lon,lat = np.meshgrid(lons,lats)
x, y = mp(lon,lat)
c_scheme =mp.pcolor(x,y , np.squeeze(eto[0,:,:]), cmap = 'jet')
mp.drawcoastlines()
mp.drawstates()
mp.drawcounties()

parallels = np.arange(30,50,5.) # make latitude lines ever 5 degrees from 30N-50N
meridians = np.arange(-95,-70,5.) # make longitude lines every 5 degrees from 95W to 70W
mp.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
mp.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

plt.title('eto')
plt.xlabel('lon')
plt.ylabel('lat')
plt.legend()
plt.grid(True)
plt.show()

data.close()
data_1.close()
data_2.close()
data_3.close()
data_4.close()
data_5.close()
data_6.close()