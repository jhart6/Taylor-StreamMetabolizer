#Mouth Metabolism Models
#All available dates
#working to eliminate as many errors as possible

library(streamMetabolizer)
library(dplyr)
library(lubridate)
library(zoo)


#####PREPARING DRIVER DATA#####
#Import data
setwd("~/Dropbox/BUTMAN/SUCR/TS Sensor Data/Metab Datasets/")
raw <- read.csv("Mouth_Metab.csv")
raw <- raw %>%
  select(DateTime.adj,Temp.C,DO.mgL,Intensity.lumft2)
long <- -122.24631

#Convert DateTime PST to Posix to Solar Time
local.time.pdt <- as.POSIXct(raw$DateTime.adj,format = "%Y-%m-%d %H:%M:%S",tz = 'America/Los_Angeles')
solar.time <- posix.solar.time <- calc_solar_time(local.time.pdt,longitude = long)
raw$DateTime.adj <- solar.time
raw$DateTime.adj <- round_date(raw$DateTime.adj,unit = "10 min")
colnames(raw) <- c("solar.time","Temp.C","DO.mgL","Intensity.lumft2")

dat <- raw

# #create vector of all possible dates
# z <- calc_solar_time(seq.POSIXt(as.POSIXct("2019-06-29"), as.POSIXct("2019-01-03"),by = "10 mins"),longitude=long)
# z <- round_date(z, unit = "10 min")
# z <- z[which(z >= as.POSIXct("2019-06-29 06:40:00") & z < as.POSIXct("2019-09-02 10:00:00"))]

#match data in 3 day dataset to derived solartime sequence
# dat <- data.frame(solar.time = z)
# dat <- left_join(dat, raw, by = "solar.time")
# 
# dat <- dat[1:9372,]

#interpolate any missing data points
# #do.obs
# plot(dat$solar.time,dat$DO.mgL)
# dat$DO.mgL <- na.approx(dat$DO.mgL)
# plot(dat$solar.time,dat$DO.mgL)

#water temp
# plot(dat$solar.time,dat$Temp.C)
# dat$Temp.C <- na.approx(dat$Temp.C)
# plot(dat$solar.time,dat$Temp.C)

#light
# dat$Intensity.lumft2 <- na.approx(dat$Intensity.lumft2)
# plot(dat$solar.time,dat$Intensity.lumft2)


#create additional necessary vars for streammetabolizer
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

#subset to start at 4AM and end at 3:50AM
# dat.mouth <- dat.mouth[which(dat.mouth$solar.time >= "2019-06-29 21:00:00" & dat.mouth$solar.time < "2019-09-01 21:00:00"),]



# #BORING MLE MODEL
# boring.name <- mm_name(type="mle")
# boring.specs <- specs(boring.name)
# boring.mm <- metab(boring.specs,dat.mouth)
# 
# get_params(boring.mm)
# plot_metab_preds(boring.mm)
# plot_DO_preds(boring.mm)
# metab.adjP <- predict_metab(boring.mm)
# do.adjP <- predict_DO(boring.mm)
# 
# param.warnings<-get_params(boring.mm) %>%
#   select(date,warnings,errors) 
# metab.warnings <- predict_metab(boring.mm) %>%
#   select(date,warnings,errors)


#BORING RESP MODEL
resp.boring.name <- mm_name(type="mle",ER_fun="q10temp")
resp.boring.specs <- specs(resp.boring.name)
resp.boring.mm <- metab(resp.boring.specs,dat.mouth)

modelparams <- get_params(resp.boring.mm)
plot_metab_preds(resp.boring.mm)
plot_DO_preds(resp.boring.mm)
metab.adjP <- predict_metab(resp.boring.mm)
do.adjP <- predict_DO(resp.boring.mm)

param.warnings<-get_params(resp.boring.mm) %>%
  select(date,warnings,errors) 
metab.warnings <- predict_metab(resp.boring.mm) %>%
  select(date,warnings,errors)


#BORING GPP MODEL
# gpp.boring.name <- mm_name(type="mle",GPP_fun = "satlight")
# gpp.boring.specs <- specs(gpp.boring.name)
# gpp.boring.mm <- metab(gpp.boring.specs,dat.mouth)
# 
# get_params(gpp.boring.mm)
# plot_metab_preds(gpp.boring.mm)
# plot_DO_preds(gpp.boring.mm)
# metab.adjP <- predict_metab(gpp.boring.mm)
# do.adjP <- predict_DO(gpp.boring.mm)
# 
# param.warnings<-get_params(gpp.boring.mm) %>%
#   select(date,warnings,errors) 
# metab.warnings <- predict_metab(gpp.boring.mm) %>%
#   select(date,warnings,errors)
# 
# 
# #GPP MODEL with AdjP+A
# gpp.name <- mm_name(type="mle",GPP_fun = "satlight")
# gpp.specs <- specs(gpp.name,init.Pmax = 8, init.alpha = 0.08)
# gpp.mm <- metab(gpp.specs,dat.mouth)
# 
# get_params(gpp.mm)
# plot_metab_preds(gpp.mm)
# plot_DO_preds(gpp.mm)
# metab.adjP <- predict_metab(gpp.mm)
# do.adjP <- predict_DO(gpp.mm)
# 
# param.warnings<-get_params(gpp.mm) %>%
#   select(date,warnings,errors) 
# metab.warnings <- predict_metab(gpp.mm) %>%
#   select(date,warnings,errors)
# 
# #SATLIGHT and Q10TEMP MODEL
# full.name <- mm_name(type="mle",GPP_fun = "satlight",ER_fun='q10temp')
# full.specs <- specs(full.name)
# full.mm <- metab(full.specs,dat.mouth)
# 
# get_params(gpp.mm)
# plot_metab_preds(gpp.mm)
# plot_DO_preds(gpp.mm)
# metab.adjP <- predict_metab(gpp.mm)
# do.adjP <- predict_DO(gpp.mm)
# 
# param.warnings<-get_params(gpp.mm) %>%
#   select(date,warnings,errors) 
# metab.warnings <- predict_metab(gpp.mm) %>%
#   select(date,warnings,errors)



#####################################################################################################
#dig a little deeper into most successful MLE model (Boring RESP)
#10/24/19

write.csv(modelparams,"boringRespParams.csv",row.names = FALSE)
setwd("~/Dropbox/BUTMAN/SUCR/StreamMetabolizer/191024/")
modelparams2 <- read.csv("boringRespParams.csv") #created up and down sd's in excel (way easier) for plotting sd of K600 estimates


#digging into K600 estimates
plot(as.Date(modelparams2$date),modelparams2$K600.daily,pch=16,ylab = "K600",xlab = "Date")
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$K600.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$K600.sd.up,code = 3,length =0)

plot(as.Date(modelparams2$date),modelparams2$GPP.daily,pch=16,ylab = "GPP",xlab = "Date",ylim=c(-10,10))
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$GPP.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$GPP.sd.up,code = 3,length =0)

plot(as.Date(modelparams2$date),modelparams2$ER20,pch=16,ylab = "ER",xlab = "Date",ylim=c(-175,5))
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$ER.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$ER.sd.up,code = 3,length =0)

quartz()
par(mfrow = c(3,1))
plot(as.Date(modelparams2$date),modelparams2$K600.daily,pch=16,ylab = "K600",xlab = "Date")
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$K600.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$K600.sd.up,code = 3,length =0)
plot(as.Date(modelparams2$date),modelparams2$GPP.daily,pch=16,ylab = "GPP",xlab = "Date",ylim=c(-10,10))
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$GPP.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$GPP.sd.up,code = 3,length =0)
plot(as.Date(modelparams2$date),modelparams2$ER20,pch=16,ylab = "ER",xlab = "Date",ylim=c(-175,5))
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$ER.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$ER.sd.up,code = 3,length =0)


#adding [CO2] dissolved to the plots
#now added [co2] to modelparams2 csv
modelparams2 <- read.csv("boringRespParams.csv") #now includes dissolved CO2 data

plot(as.Date(modelparams2$date),modelparams2$ER20,pch=16,ylab = "ER",xlab = "Date",ylim=c(-175,5))
abline(0,0,lty=2)
arrows(x0 = as.Date(modelparams2$date), y0 = modelparams2$ER.sd.down, x1 = as.Date(modelparams2$date), y1 = modelparams2$ER.sd.up,code = 3,length =0)
points(as.Date(modelparams2$date),modelparams2$X)


