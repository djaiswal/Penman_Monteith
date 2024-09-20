

ET <- function(tmean, wind_speed, 
               relative_humidity, net_radiation, 
               P =101325,data_co2,co2_ref = 330,G=0,gs_ref=0.0061 )
  {
  # Calculate daily mean temperature (Tmean) in Celsius
  # tmean = ((temperature_max + temperature_min) / 2.0)
  #net_radiation = (long_radiation + short_radiation)*(0.77)
  gs = gs_ref * (1 / (1 + 0.663 * ((data_co2 / co2_ref) - 1)))
  Rn = net_radiation
  
  # Calculate saturation vapor pressure (es) and actual vapor pressure (ea)
  es = 0.6108 * exp((17.27 * tmean) / (tmean + 237.3))
  ea = (relative_humidity / 100.0) * es
  
  # Calculate the slope of the saturation vapor pressure curve (delta)
  delta = (4098 * es) / ((tmean + 237.3) ** 2)
  
  # Calculate the psychrometric constant (gamma)
  gamma = 0.000655*P
  
  # Calculate ET0 using the Penman-Monteith equation
  et0 = ((0.408 * delta * (Rn - G) + 
            (gamma * (900 / (tmean + 273)) * wind_speed * (es - ea))) /       
           (delta + gamma * (1 + 0.035 * wind_speed)/gs))
}



value <- ET(tmean = 32, wind_speed=1.36, 
            relative_humidity=57.1, net_radiation = 41.22, 
            P =77325,data_co2=1012,co2_ref = 330)
value

#Temp Sensitivity
input_data_Temp_sensitivity <- data.frame(tmean=seq(-10,50,length =20),
                                          wind_speed = rep(1.36,20),
                                          relative_humidity=rep(57.1,20),
                                          net_radiation=rep(53,20),
                                          data_co2 = rep(1012,20))

ETout <- numeric(20)
for(i in 1:20){
  ETout[i] = ET(tmean = input_data_Temp_sensitivity$tmean[i], 
                   wind_speed=input_data_Temp_sensitivity$wind_speed[i], 
                   relative_humidity=input_data_Temp_sensitivity$relative_humidity[i], 
                   net_radiation = input_data_Temp_sensitivity$net_radiation[i], 
                   P =101325,
                   data_co2=input_data_Temp_sensitivity$data_co2[i],
                   co2_ref = 330)
}

Tem_Sensitivity_data <- data.frame(Temp = input_data_Temp_sensitivity$tmean, ET=ETout)
xyplot(ET~Temp, data=Tem_Sensitivity)

#Rad Sensitivity
input_data_Rad_sensitivity <- data.frame(tmean=rep(30,20),
                                          wind_speed = rep(1.36,20),
                                          relative_humidity=rep(57.1,20),
                                          net_radiation=seq(20,60,length=20),
                                          data_co2 = rep(1012,20))
ETout <- numeric(20)
for(i in 1:20){
  ETout[i] = ET(tmean = input_data_Rad_sensitivity$tmean[i], 
                     wind_speed=input_data_Rad_sensitivity$wind_speed[i], 
                     relative_humidity=input_data_Rad_sensitivity$relative_humidity[i], 
                     net_radiation = input_data_Rad_sensitivity$net_radiation[i], 
                     P =101325,
                     data_co2=input_data_Rad_sensitivity$data_co2[i],
                     co2_ref = 330)
}

Rad_Sensitivity_data <- data.frame(Rad = input_data_Rad_sensitivity$net_radiation, ET=ETout)
xyplot(ET~Rad, data=Rad_Sensitivity_data)


#CO2 Sensitivity
input_data_Rad_sensitivity <- data.frame(tmean=rep(30,20),
                                         wind_speed = rep(1.36,20),
                                         relative_humidity=rep(57.1,20),
                                         net_radiation=rep(53,20),
                                         data_co2 = seq(400,1200,length=20))
ETout <- numeric(20)
for(i in 1:20){
  ETout[i] = ET(tmean = input_data_Rad_sensitivity$tmean[i], 
                wind_speed=input_data_Rad_sensitivity$wind_speed[i], 
                relative_humidity=input_data_Rad_sensitivity$relative_humidity[i], 
                net_radiation = input_data_Rad_sensitivity$net_radiation[i], 
                P =101325,
                data_co2=input_data_Rad_sensitivity$data_co2[i],
                co2_ref = 330)
}

CO2_Sensitivity_data <- data.frame(CO2= input_data_Rad_sensitivity$data_co2, ET=ETout)
xyplot(ET~CO2, data=CO2_Sensitivity_data)

#Windspeed Sensitivity
input_data_ws_sensitivity <- data.frame(tmean=rep(30,20),
                                         wind_speed = seq(0.4,4,length=20),
                                         relative_humidity=rep(57.1,20),
                                         net_radiation=rep(53,20),
                                         data_co2 = rep(1012,20))
ETout <- numeric(20)
for(i in 1:20){
  ETout[i] = ET(tmean = input_data_ws_sensitivity$tmean[i], 
                wind_speed=input_data_ws_sensitivity$wind_speed[i], 
                relative_humidity=input_data_ws_sensitivity$relative_humidity[i], 
                net_radiation = input_data_ws_sensitivity$net_radiation[i], 
                P =101325,
                data_co2=input_data_ws_sensitivity$data_co2[i],
                co2_ref = 330)
}

ws_Sensitivity_data <- data.frame(ws= input_data_ws_sensitivity$wind_speed, ET=ETout)
xyplot(ET~ws, data=ws_Sensitivity_data)

#RH Sensitivity
input_data_rh_sensitivity <- data.frame(tmean=rep(30,20),
                                         wind_speed = rep(1.36,20),
                                         relative_humidity=seq(26,96,length=20),
                                         net_radiation=rep(53,20),
                                         data_co2 = rep(1012,20))
ETout <- numeric(20)
for(i in 1:20){
  ETout[i] = ET(tmean = input_data_rh_sensitivity$tmean[i], 
                wind_speed=input_data_rh_sensitivity$wind_speed[i], 
                relative_humidity=input_data_rh_sensitivity$relative_humidity[i], 
                net_radiation = input_data_rh_sensitivity$net_radiation[i], 
                P =101325,
                data_co2=input_data_rh_sensitivity$data_co2[i],
                co2_ref = 330)
}

rh_Sensitivity_data <- data.frame(rh= input_data_rh_sensitivity$relative_humidity, ET=ETout)
xyplot(ET~rh, data=rh_Sensitivity_data)

