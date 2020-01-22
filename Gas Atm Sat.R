###script to estimate gas saturation given equilibrium with atmosphere
##21 January 2020
##Adapted from Luke Loken
#JAH

#1) need Henry's constant (gas-specific)
#can do this two ways: Plummer and Busenberg 1982; Geochimica et Cosmochimica Acta
#2) get atmospheric saturation based on Henry's constant, atm pressure, and gas (gas-specific)

#function to calculate Henry's constant from temperature for CO2 (Plummer and Busenberg)
#temperature needs to be in Kelvin
Kh_Plummer <- function(Temperature){
  kh_Plummer=10^(108.3865 + 0.01985076*Temperature - 6919.53*(Temperature^-1) -40.45154*log10(Temperature)+669365*(Temperature^-2))
  return(kh_Plummer)
}    

#function to calculate Henry's constant from temperature for many gases
getKh <- function(temperature, gas){
  Kh  <-  data.frame("O2"=c(1.3*10^-3, 1700))
  Kh <- cbind(Kh,"H2"=c(7.8*10^-4,500))
  Kh <- cbind(Kh,"CO2"= c(3.4*10^-2, 2400 ))
  Kh <- cbind(Kh,"N2"=c(6.1*10^-4, 1300))
  Kh <- cbind(Kh,"He"=c(3.7*10^-4, 230))
  Kh <- cbind(Kh,"Ne"=c(4.5*10^-4,490))
  Kh <- cbind(Kh,"Ar"=c(1.4*10^-3, 1300))
  Kh <- cbind(Kh,"CO"=c(9.5*10^-4,1300))
  Kh <- cbind(Kh, "O3"=c(1.2*10^-2, 2300))
  Kh <- cbind(Kh, "N2O"=c(2.5*10^-2, 2600))
  Kh <- cbind(Kh, "SF6"=c(2.4*10^-4, 2400))
  Kh <- cbind(Kh, "CH4"=c(1.4*10^-3, 1700))
  Kh <- cbind(Kh, "C3H8"=c(1.4*10^-3, 2700))
  Kh <- cbind(Kh, "NH3"=c(5.6*10^1, 4200))
  
  if (!is.character(gas)){stop(paste('gas must be a character. was given as',gas))}
  if (!any(names(Kh)==gas)){stop(paste(gas,'not found in list of coded gases'))}
  
  Khprime <-  unlist(Kh[gas])[1]
  C  <-  unlist(Kh[gas])[2]
  
  LakeKh= as.numeric(Khprime*exp(C*((1/temperature)-(1/298))))
}

#function to get atmospheric saturation based on Henry's constant, atm pressure, and gas
getSaturation <- function(LakeKh, AtmP, gas){
  Atmosphere  <-  data.frame("O2"=209460)
  Atmosphere <- cbind(Atmosphere,"H2"=0.55)
  Atmosphere <- cbind(Atmosphere, "N2"=780840)
  Atmosphere <- cbind(Atmosphere, "Ar"=9340)
  Atmosphere <-cbind(Atmosphere, "CO2"=400)
  Atmosphere <-cbind(Atmosphere, "He"=5.24)
  Atmosphere <-cbind(Atmosphere, "Ne"=18.18)
  Atmosphere <- cbind(Atmosphere, "CH4"= 1.91)
  Atmosphere <- cbind(Atmosphere, "O3"=0.07)#potential max concentration
  Atmosphere <-cbind(Atmosphere, "N2O"= 0.325)
  Atmosphere <-cbind(Atmosphere, "CO"=0.1)
  Atmosphere <-cbind(Atmosphere, "NH3"=NA)
  
  if (!is.character(gas)){stop(paste('gas must be a character. was given as',gas))}
  if (!any(names(Atmosphere)==gas)){stop(paste(gas,'not found in list of coded gases'))}
  
  AtmosphereConc <-  unlist(Atmosphere[gas])[1]
  
  EquilSaturation=AtmosphereConc*LakeKh/AtmP #umol/L, mmol/m3
}
#produces equilibrium saturation in uM (micromolar)


Elevation = 16
Pressure=(1-(.0000225577*Elevation))^5.25588

