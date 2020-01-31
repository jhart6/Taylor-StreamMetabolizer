#Getting to Know StreamMetabolizer
#Working with time series sensor data from Mouth site on Taylor Creek
#19 September 2019
#JAH

install.packages("streamMetabolizer", dependencies=TRUE, 
                 repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

library(streamMetabolizer)
library(dplyr)
library(lubridate)
library(zoo)

#####PREPARING DRIVER DATA#####
#Import mouth time series dataset
setwd("~/Dropbox/BUTMAN/SUCR/TS Sensor Data/Metab Datasets/")
raw <- read.csv("Mouth_Metab.csv")
raw <- raw %>% #select only desired columns
  select(DateTime.adj,Temp.C,DO.mgL,Intensity.lumft2)
long <- -122.24631 #longitude of field site

#Convert DateTime PST to Posix to Solar Time (package requires Solar Time)
local.time.pdt <- as.POSIXct(raw$DateTime.adj,format = "%Y-%m-%d %H:%M:%S",tz = 'America/Los_Angeles')
solar.time <- posix.solar.time <- calc_solar_time(local.time.pdt,longitude = long)

#Estimate DO.sat (DO concentration at equilibrium saturation)
DO.sat <- calc_DO_sat(temp = raw$Temp.C,press = 1013, sal = 0)

#Convert Intensity.lumft2 to correct light units (umol/m^2/s)
light <- raw$Intensity.lumft2 * 10.7639104167 * 0.0185

#Create a depth column (m)
depth = 0.02 #value estimated based on sensor placement in Taylor Creek

#Create dataframe with relevant Metab data in correct units
dat <- data.frame(solar.time,DO.obs = raw$DO.mgL,DO.sat,temp.water = raw$Temp.C,light,depth = rep(depth,length(solar.time)))

#subset data to only dates with available light data
#starts 6/28/19 at 12:04 PST
dat <- dat[which(is.na(dat$light)==F),]



#####INSPECTING DRIVER DATA#####
#plot all three time series
plot(dat$solar.time,dat$DO.obs)
plot(dat$solar.time,dat$temp.water,type='l')
plot(dat$solar.time,dat$light,type='l')

#three-panel plot with same x-axis
quartz()
par(mfrow=c(3,1),mar=c(0,4,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
plot(dat$solar.time,dat$DO.obs,xlab=NA,type='l',xaxt='n',ylab="DO (mg/l)",col='red')

par(mar=c(0,4,0,4),mgp=c(1.5,0.5,0),tck=-0.02)
plot(dat$solar.time,dat$temp.water,type='l',xlab=NA,xaxt='n',ylab="Water Temp (C)",col='blue')

par(mar=c(3,4,0,4),mgp=c(1.5,0.5,0),tck=-0.02)
plot(dat$solar.time,dat$light,type='l',ylab="Light (umol/m^2/s)", xlab="Date",col='black')


# #cross reference with met data
###1/31/20 - seatac_met.csv is not currently the same length as mouth_metab
# met <- read.csv('seatac_met.csv')
# 
# plot(as.Date(met$DATE),met$PRCP,type='l')
# plot(as.Date(met$DATE),met$TAVG,type='l')
# 
# 
# #put met and stream driver data together in the same plot
# quartz()
# par(mfrow=c(5,1),mar=c(0,4,1,4),mgp=c(1.5,0.5,0),tck=-0.02,bg='white')
# plot(as.Date(met$DATE),met$PRCP,type='l',xaxt='n',ylab='Precipitation (in)')
# par(mar=c(0,4,0,4),mgp=c(1.5,0.5,0),tck=-0.02)
# plot(as.Date(met$DATE),met$TAVG,type='l',ylab='Average Air Temp (F)',xlab="Date")
# 
# par(mar=c(0,4,3,4),mgp=c(1.5,0.5,0),tck=-0.02)
# plot(dat$solar.time,dat$DO.obs,xlab=NA,type='l',xaxt='n',ylab="DO (mg/l)",col='red')
# par(mar=c(0,4,0,4),mgp=c(1.5,0.5,0),tck=-0.02)
# plot(dat$solar.time,dat$temp.water,type='l',xlab=NA,xaxt='n',ylab="Water Temp (C)",col='blue')
# par(mar=c(3,4,0,4),mgp=c(1.5,0.5,0),tck=-0.02)
# plot(dat$solar.time,dat$light,type='l',ylab="Light (umol/m^2/s)", xlab="Date",col='black')


#####The following model setups are varied according to GPP and ER functions. Different combinations result in different metabolism estimates, with various levels of success.


#####DEFAULT METAB MODEL#####
#model = mle1 (basic, default MLE metab model)
mle1.name <- mm_name(type = "mle") #define model structure (basic mle model, default settings)
mle1.specs <- specs(mle1.name) #assign an object to the model specs for inspection
mm1 <- metab(mle1.specs, data = dat) #run metabolism model

#save this model's output
daily.metab.mle1 <- predict_metab(mm1)
daily.do.mle1 <- predict_DO(mm1)

#plot this model's output
quartz()
plot_metab_preds(mm1)
quartz()
plot_DO_preds(mm1)




####SATURATOR METAB MODEL####
#model = mle.gppfun.satlight
mle2.name = mm_name(type = "mle", GPP_fun = 'satlight',ER_fun = 'constant')
mle2.specs = specs(mle2.name)
mm2 <- metab(mle2.specs,data = dat)

#save this model's output
daily.metab.mle2 <- predict_metab(mm2)
daily.do.mle2 <- predict_DO(mm2)

#plot this model's output
quartz()
plot_metab_preds(mm2)
quartz()
plot_DO_preds(mm2)

#look at model prediction warning messages
get_params(mm2) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params

predict_metab(mm2) %>%
  select(date,warnings,errors)




#####SATURATOR 2.0#####
#use data_daily to make Pmax and alpha estimates from the previous model (Saturator)
#dynamic Pmax and alpha estimates (rather than fixed)

#model = mle.gppfun.satlight
mle3.name = mm_name(type = "mle", GPP_fun = 'satlight',ER_fun = 'constant')
mle3.specs = specs(mle3.name)
mm3 <- metab(mle3.specs,data = dat,data_daily = select(get_params(mm2),date,init.Pmax = Pmax, init.alpha = alpha))

#save this model's output
daily.metab.mle3 <- predict_metab(mm3)
daily.do.mle3 <- predict_DO(mm3)

#plot this model's output
quartz()
plot_metab_preds(mm3)
quartz()
plot_DO_preds(mm3)

#look at model prediction warning messages
get_params(mm3) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params

predict_metab(mm3) %>%
  select(date,warnings,errors) #warnings about metabolism estimates




#####SATURATOR 3.0#####
#update Pmax and alpha started values based on previous model (Saturator)
#fixed, but updated Pmax and alpha estimates
#Pmax = 6.4
#alpha = 0.008

#model = mle.gppfun.satlight
mle4.name = mm_name(type = "mle", GPP_fun = 'satlight',ER_fun = 'constant')
mle4.specs = specs(mle4.name, init.Pmax = 6.4, init.alpha = 0.008)
mm4 <- metab(mle4.specs,data = dat)

#save this model's output
daily.metab.mle4 <- predict_metab(mm4)
daily.do.mle4 <- predict_DO(mm4)

#plot this model's output
quartz()
plot_metab_preds(mm4)
quartz()
plot_DO_preds(mm4)

#look at model prediction warning messages
get_params(mm4) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params

predict_metab(mm4) %>%
  select(date,warnings,errors) #warnings about metabolism estimates




#####SATURATOR 4.0#####
#MIXED date-specific and fixed Pmax and alpha values
#fixed values are the default values

#model = mle.gppfun.satlight
mle5.name = mm_name(type = "mle", GPP_fun = 'satlight',ER_fun = 'constant')
mle5.specs = specs(mle5.name, init.Pmax = 10, init.alpha = 0.0001)
mm5 <- metab(mle5.specs,data = dat)

#save this model's output
daily.metab.mle5 <- predict_metab(mm5)
daily.do.mle5 <- predict_DO(mm5)

#plot this model's output
quartz()
plot_metab_preds(mm5)
quartz()
plot_DO_preds(mm5)

#look at model prediction warning messages
get_params(mm5) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params

predict_metab(mm4) %>%
  select(date,warnings,errors) #warnings about metabolism estimates



#####Saturator + Q10 Resp#####
mle6.name = mm_name(type = "mle", GPP_fun = 'satlight',ER_fun = 'q10temp')
mle6.specs = specs(mle6.name)
mm6 <- metab(mle6.specs,data = dat)

#save this model's output
daily.metab.mle6 <- predict_metab(mm6)
daily.do.mle6 <- predict_DO(mm6)

#plot this model's output
quartz()
plot_metab_preds(mm6)
quartz()
plot_DO_preds(mm6)

#look at model prediction warning messages
get_params(mm6) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params

predict_metab(mm6) %>%
  select(date,warnings,errors) #warnings about metabolism estimates




#####Saturator + Q10 Resp with Dynamic Init Values#####
#use data_daily to make Pmax and alpha estimates from the previous model (Saturator + Q10)
#dynamic Pmax and alpha estimates (rather than fixed)

mle7.name = mm_name(type = "mle", GPP_fun = 'satlight',ER_fun = 'q10temp')
mle7.specs = specs(mle7.name)
mm7 <- metab(mle7.specs,data = dat,data_daily = select(get_params(mm6),date,init.Pmax = Pmax, init.alpha = alpha))

#save this model's output
daily.metab.mle7 <- predict_metab(mm7)
daily.do.mle7 <- predict_DO(mm7)

#plot this model's output
quartz()
plot_metab_preds(mm7)
quartz()
plot_DO_preds(mm7)

#look at model prediction warning messages
get_params(mm7) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params

predict_metab(mm7) %>%
  select(date,warnings,errors) #warnings about metabolism estimates


#####Create P-I Plot#####
#calculate daily mean surface irradiance
mean_irr <- dat %>%
  select(solar.time,light) %>%
  group_by(as.Date(solar.time))%>%
  summarise(avg = mean(light)) 
colnames(mean_irr) <- c("date","mean.irr")

gpp.satlight.q10 <- daily.metab.mle6$GPP
gpp.satlight.constant <- daily.metab.mle2$GPP  
gpp.linlight.constant <- daily.metab.mle1$GPP

mod1 <- lm(gpp.satlight.q10~mean_irr$mean.irr)
mod2 <- lm(gpp.satlight.constant~mean_irr$mean.irr)
mod3 <- lm(gpp.linlight.constant~mean_irr$mean.irr)

quartz()
par(mfrow=c(3,1))
xlab = "Daily Mean Surface Light"
plot(mean_irr$mean.irr,gpp.satlight.q10,pch=16,xlab=xlab)
abline(mod1)
plot(mean_irr$mean.irr,gpp.satlight.constant,pch=16,xlab=xlab)
abline(mod2)
plot(mean_irr$mean.irr,gpp.linlight.constant,pch=16,xlab=xlab)
abline(mod3)




#####Get MLE to Work for 3 Days of Data####
#isolate three days of the data
dat.Aug1 <- dat[which(date(dat$solar.time)=="2019-08-01"),]
dat.Aug2 <- dat[which(date(dat$solar.time)=="2019-08-02"),]
dat.Aug3 <- dat[which(date(dat$solar.time)=="2019-08-03"),]
dat.Aug4 <- dat[which(date(dat$solar.time)=="2019-08-04"),]
dat.3day <- rbind(dat.Aug1,dat.Aug2,dat.Aug3,dat.Aug4)

#round seconds down to the minute
dat.3day$solar.time <- round_date(dat.3day$solar.time,unit = "minute")
dat.3day <- dat.3day[which(dat.3day$solar.time >= "2019-07-31 21:00:00" & dat.3day$solar.time < "2019-08-03 21:00:00"),]

#create solartime sequence for all desired time points
z <- calc_solar_time(seq.POSIXt(as.POSIXct("2019-08-01"), as.POSIXct("2019-08-05"),by = "10 min"),longitude=long)
z <- round_date(z, unit = "10 min")
z <- z[which(z >= as.POSIXct("2019-07-31 21:00:00") & z < as.POSIXct("2019-08-03 21:00:00"))]

#match data in 3 day dataset to derived solartime sequence
three.dat <- data.frame(solar.time = z)
three.dat <- left_join(three.dat, dat.3day, by = "solar.time")

#interpolate any missing data points
#do.obs
plot(three.dat$solar.time,three.dat$DO.obs)
three.dat$DO.obs <- na.approx(three.dat$DO.obs)
plot(three.dat$solar.time,three.dat$DO.obs)

#water temp
plot(three.dat$solar.time,three.dat$temp.water)
three.dat$temp.water <- na.approx(three.dat$temp.water)
plot(three.dat$solar.time,three.dat$temp.water)

#do.sat
three.dat$DO.sat <- calc_DO_sat(temp = three.dat$temp.water,press = 1013, sal = 0)
plot(three.dat$solar.time,three.dat$DO.sat)

#light
three.dat$light <- na.approx(three.dat$light)
plot(three.dat$solar.time,three.dat$light)

#depth
three.dat$depth = rep(0.2, length(three.dat$solar.time))

#run a basic mle metabolism model (satlight, q10temp) on 3 days
mle.3day.name <- mm_name(type="mle",GPP_fun = "satlight", ER_fun = "q10temp")
mle.3day.specs <- specs(mle.3day.name,init.Pmax = 7.2, init.alpha = 0.01)
mm.mle.3day <- metab(mle.3day.specs,three.dat)

plot_metab_preds(mm.mle.3day)
plot_DO_preds(mm.mle.3day)
metab.3day <- predict_metab(mm.mle.3day)
do.3day <- predict_DO(mm.mle.3day)

param.warnings<-get_params(mm.mle.3day) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(mm.mle.3day) %>%
  select(date,warnings,errors) 




#####SIMULATION#####
sim.name <- mm_name("sim",GPP_fun = "satlight", ER_fun = "q10temp")
mm.sim.trial <- metab(specs(sim.name),three.dat)
plot_DO_preds(mm.sim.trial)
plot_metab_preds(mm.sim.trial)

params.sim.trial <- data.frame(date = as.Date(paste0("2019-08-",01:03)),GPP.daily = 2.1, ER20 = -5:-3, K600.daily=16)
mm.sim.trial2 <- metab(specs(sim.name),three.dat,data_daily = params.sim.trial)
plot_DO_preds(mm.sim.trial2)
plot_metab_preds(mm.sim.trial2)

params.sim.prev <- get_params(mm.mle.3day, uncertainty = "none", messages = FALSE)
mm.sim.prev <- metab(specs(sim.name),three.dat,data_daily = params.sim.prev)
plot_DO_preds(mm.sim.prev)
plot_metab_preds(mm.sim.prev)

