#Script to explore oxygen dynamics at Mouth



##################################
####Load Necessary Packages#######
##################################
library(streamMetabolizer)
library(dplyr)
library(lubridate)
library(zoo)


##################################
###########Import data############
##################################
setwd("~/Dropbox/GitHub Repos/Taylor-StreamMetabolizer/data/Metab Datasets/")
raw <- read.csv("Mouth_Metab.csv")
raw <- raw %>%
  select(DateTime.adj,Temp.C,DO.mgL,Intensity.lumft2)
long <- -122.24631


#####################################################
######Convert DateTimePST to Posix to Solar Time#####
#####################################################
local.time.pdt <- as.POSIXct(raw$DateTime.adj,format = "%Y-%m-%d %H:%M:%S",tz = 'America/Los_Angeles')
solar.time <- posix.solar.time <- calc_solar_time(local.time.pdt,longitude = long)
raw$DateTime.adj <- solar.time
raw$DateTime.adj <- round_date(raw$DateTime.adj,unit = "10 min")
colnames(raw) <- c("solar.time","Temp.C","DO.mgL","Intensity.lumft2")

dat <- raw


##########################################
#####create additional necessary vars#####
##########################################
#depth
depth <- rep(0.2,length(dat$solar.time))

#Estimate DO.eq (DO concentration at equilibrium saturation in mg/L)
DO.eq <- calc_DO_sat(temp = dat$Temp.C,press = 1013, sal = 0)

#Estimate DO.sat (DO % saturation)
DO.sat <- (dat$DO.mgL/DO.eq) * 100

#Convert Intensity.lumft2 to correct light units (umol/m^2/s)
light <- dat$Intensity.lumft2 * 10.7639104167 * 0.0185

#Create dataframe with relevant Metab data in correct units
dat.mouth <- data.frame(solar.time=dat$solar.time,DO.obs = dat$DO.mgL,DO.sat,temp.water = dat$Temp.C,light,depth = depth)


#########################################
######Plot DO mg.L and DO sat############
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(dat.mouth$solar.time,dat.mouth$DO.obs,pch=1,col="blue",xlab="Date",ylab="DO.obs (mg/L)")
par(new=T)
plot(dat.mouth$solar.time,dat.mouth$DO.sat,axes=F,xlab=NA,ylab=NA,pch=1,col="red")
axis(side=4)
mtext(side=4,line=2,"DO %Sat")


#########################################
######Plot DO % sat######################
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(dat.mouth$solar.time,dat.mouth$DO.sat,col="red",xlab="Date",ylab="DO %Sat")
abline(100,0,col="black",lty=2,lwd=2)


#########################################
######Load in EXO2 Oxygen Data##########
#########################################
setwd("~/Dropbox/GitHub Repos/Taylor-StreamMetabolizer/data/")
exo<-read.csv("EXO2 Oxygen Data.csv")

local.time.pdt.exo <- as.POSIXct(exo$DateTime.Rounded,format = "%Y-%m-%d %H:%M:%S",tz = 'America/Los_Angeles')
solar.time.exo <- posix.solar.time.exo <- calc_solar_time(local.time.pdt.exo,longitude = long)
exo$solar.time <- solar.time.exo
exo$solar.time <- round_date(exo$solar.time,unit = "10 min")


#########################################
######Compare DO mg.L####################
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(dat.mouth$solar.time,dat.mouth$DO.obs,pch=1,col="blue",xlab="Date",ylab="DO.mgL MiniDOT",ylim=c(0,14))
par(new=T)
plot(as.Date(exo$solar.time),exo$odo.mgL,pch=16,cex=2,axes=F,ylab=NA,xlab=NA,ylim=c(0,14))
axis(side=4)
mtext(side=4,line=2,"DO.mgL EXO2")
legend("bottomright",c("MiniDOT TS mg.L","EXO2 mg.L"),pch=c(1,16),col=c("blue","black"))


#########################################
######Compare DO sat####################
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(dat.mouth$solar.time,dat.mouth$DO.sat,pch=1,col="red",xlab="Date",ylab="DO.sat MiniDOT",ylim=c(0,120))
par(new=T)
plot(as.Date(exo$solar.time),exo$odo.sat,pch=16,cex=2,axes=F,ylab=NA,xlab=NA,ylim=c(0,120))
axis(side=4)
mtext(side=4,line=2,"DO.sat EXO2")
legend("bottomright",c("MiniDOT TS DO.sat","EXO2 DO.sat"),pch=c(1,16),col=c("red","black"))


#########################################
#####Create Dataframe to Compare Numbers#
#########################################
dat.mouth.subset = dat.mouth %>%
  select(solar.time,DO.obs,DO.sat)

exo.mouth.subset = exo %>%
  select(solar.time,odo.mgL,odo.sat)

combo = left_join(exo.mouth.subset,dat.mouth.subset,by="solar.time")


#########################################
#####Estimate Correction for mg.L########
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(combo$odo.mgL,combo$DO.obs,xlim=c(7.5,13),ylim=c(7.5,13),pch=16,xlab="DO.mgL EXO2",ylab="DO.mgL MiniDOT")
mgL.mod<-lm(combo$DO.obs~combo$odo.mgL)
abline(0,1,col="red")

quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
residual.mgL <- combo$DO.obs - combo$odo.mgL
plot(combo$odo.mgL,residual.mgL,pch=16,xlab="DO.mgL EXO2",ylab="Residual")
abline(0,0,lty=2,col="red")
residmod.mgL <- lm(residual.mgL~combo$odo.mgL)
abline(residmod.mgL)

#########################################
#####Estimate Correction for sat########
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(combo$odo.sat,combo$DO.sat,xlim=c(0,120),ylim=c(0,120),pch=16,xlab="DO.sat EXO2",ylab="DO.sat MiniDOT")
sat.mod<-lm(combo$DO.sat~combo$odo.sat)
abline(0,1,col="red")

quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
residual.sat <- combo$DO.sat - combo$odo.sat
plot(combo$odo.sat,residual.sat,pch=16,xlab="DO.sat EXO2",ylab="Residual")
abline(0,0,lty=2,col="red")
residmod.sat <- lm(residual.sat~combo$odo.sat)
abline(residmod.sat)



#########################################
#####Start adding CO2 data for comps#####
#########################################
setwd("~/Dropbox/GitHub Repos/Taylor-StreamMetabolizer/data/")
co2 <- read.csv("CO2 data.csv")
temp.obs <- co2$Temp.C

local.time.pdt.co2 <- as.POSIXct(co2$DateTime.Rounded,format = "%Y-%m-%d %H:%M:%S",tz = 'America/Los_Angeles')
solar.time.co2 <- posix.solar.time.co2 <- calc_solar_time(local.time.pdt.co2,longitude = long)
co2$solar.time <- solar.time.co2
co2$solar.time <- round_date(co2$solar.time,unit = "10 min")

#create vector of CO2 obs in ppm
co2.obs <- co2$CO2.ppm 

#create vector of CO2 sat in ppm
source("~/Dropbox/GitHub Repos/Taylor-StreamMetabolizer/Gas Atm Sat.R")
obs_henry_co2 <- Kh_Plummer(temp.obs + 273.15)
obs_CO2sat_uM <- getSaturation(LakeKh = obs_henry_co2, AtmP = Pressure, gas = "CO2")
obs_CO2sat_ppm <- obs_CO2sat_uM/obs_henry_co2 * Pressure

#estimate CO2 saturation
co2.sat <- (co2.obs/obs_CO2sat_ppm) * 100
co2 <- cbind(co2,co2.sat)

#create dataframe to match all solar.time dates
co2.mouth.subset = co2 %>%
  select(solar.time, CO2.ppm, co2.sat)

dat.mouth.subset = left_join(dat.mouth.subset, exo.mouth.subset, by = "solar.time")
dat.mouth.subset = left_join(dat.mouth.subset, co2.mouth.subset, by = "solar.time")




#########################################
#####Add CO2 Saturation Data to DO Plot##
#########################################
quartz()
par(mar=c(3,3,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(dat.mouth.subset$solar.time,dat.mouth.subset$DO.sat,pch=1,col="lavenderblush4",xlab="Date",ylab="DO.sat MiniDOT",ylim=c(0,120))
par(new=T)
plot(dat.mouth.subset$solar.time,dat.mouth.subset$odo.sat,pch=16,cex=2,axes=F,ylab=NA,xlab=NA,ylim=c(0,120),col="red")
axis(side=4)
mtext(side=4,line=2,"% Saturation")
par(new=T)
plot(dat.mouth.subset$solar.time,dat.mouth.subset$co2.sat,pch=16,cex=2,axes=F,ylab=NA,xlab=NA,ylim=c(0,120),col="blue")
legend("bottomright",c("MiniDOT TS DO % Sat","EXO2 DO % Sat","CO2 % Sat"),pch=c(1,16,16),col=c("lavenderblush4","red","blue"),cex=1.5)



