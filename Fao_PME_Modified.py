import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def fao56_penman_monteith_modified(net_rad, tmax, tmin, ws, rh_mean, atmos_pres, con_co2, shf=0.0, gs_ref=0.0061, co2_ref=330):
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
        es = round(es, 2)  # Round to two decimal points
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

# Read the Excel file(WeatherForPyeto)
data = pd.read_excel('.\WeatherForPyeto.xlsx')
# Reading The Excel Sheet Rh_Eto
Rh_Eto_data = pd.read_excel('.\WeatherForPyeto.xlsx' , sheet_name='Rh_Eto')

# Reading The Excel Sheet Tmin_Eto
Tmin_Eto_data = pd.read_excel('.\WeatherForPyeto.xlsx' , sheet_name='Tmin_Eto')

# Reading The Excel Sheet Tmax_Eto
Tmax_Eto_data = pd.read_excel('.\WeatherForPyeto.xlsx' , sheet_name='Tmax_Eto')

net_rad = data['netsolar'].values
ws = data['windspeed (m/s)'].values
temperature_min = data['Min air Temp in celsius'].values
temperature_max = data['Max temp (deg celsius)'].values
rh_mean = data['Humidity in %'].values
atmos_pres = data['atmp'].values
con_co2 = data['con(co2)ppm'].values

#Getting the data from Excel sheet (Rh_Eto)
RH = Rh_Eto_data['Humidity in %'].values

#Getting the data from Excel sheet (Tmin_Eto)
Tmin = Tmin_Eto_data['Min air Temp in celsius'].values

#Getting the data from Excel sheet (Tmax_Eto)
Tmax = Tmax_Eto_data['Max temp (deg celsius)'].values


# Example usage of fao56_penman_monteith_modified function
eto = fao56_penman_monteith_modified(net_rad, temperature_max, temperature_min, ws, rh_mean, atmos_pres, con_co2)
print("Reference Evapotranspiration (ETo):", eto)

# Create a line graph for Co2 Concentration and Evapotranspiration 
plt.plot(con_co2, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Co2 Concentration (ppm)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Co2 Concentration vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()

# Create a line graph for Relative Humidity and Evapotranspiration 
plt.plot(RH, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Relative Humidity (%)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Relative Humidity vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()

# Create a line graph for Tmin and Evapotranspiration 
plt.plot(Tmin, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Minimum Temperature (Degree Celsius)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Minimum Temperature  vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()

# Create a line graph for Tmax and Evapotranspiration 
plt.plot(Tmax, eto, label=' Rate of Evapotranspiration')
plt.xlabel('Maximum Temperature (Degree Celsius)')
plt.ylabel('Evapotranspiration (mm/day)')
plt.title('Maximum Temperature  vs Evapotranspiration')
plt.legend()
plt.grid(True)
plt.show()

