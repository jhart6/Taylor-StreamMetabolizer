#Build MLE Model for 3 days of data
  #interpolate missing data to achieve even timesteps
#Run that model on all days of data
  #interpolate missing data to achieve even timesteps


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

#Estimate DO.sat (DO concentration at equilibrium saturation)
DO.sat <- calc_DO_sat(temp = raw$Temp.C,press = 1013, sal = 0)

#Convert Intensity.lumft2 to correct light units (umol/m^2/s)
light <- raw$Intensity.lumft2 * 10.7639104167 * 0.0185

#Create a depth column (m)
depth = 0.02

#Create dataframe with relevant Metab data in correct units
dat <- data.frame(solar.time,DO.obs = raw$DO.mgL,DO.sat,temp.water = raw$Temp.C,light,depth = rep(depth,length(solar.time)))

#subset data to only dates with available light data
#starts 6/28/19 at 12:04 PST
dat <- dat[which(is.na(dat$light)==F),]





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



#####Saturator + Q10Temp MLE Model#####
#run a basic mle metabolism model (satlight, q10temp) on 3 days
default.name <- mm_name(type="mle",GPP_fun = "satlight", ER_fun = "q10temp")
default.specs <- specs(default.name)
default.mm <- metab(default.specs,three.dat)

plot_metab_preds(default.mm)
plot_DO_preds(default.mm)
metab.default <- predict_metab(default.mm)
do.default <- predict_DO(default.mm)

param.warnings<-get_params(default.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(default.mm) %>%
  select(date,warnings,errors) 

#Adjust init.Pmax until all three days work
adjP.name <- mm_name(type="mle",GPP_fun = "satlight", ER_fun = "q10temp")
adjP.specs <- specs(adjP.name,init.Pmax = 8,init.alpha = 0.08)
adjP.mm <- metab(adjP.specs,three.dat)

get_params(adjP.mm)

plot_metab_preds(adjP.mm)
plot_DO_preds(adjP.mm)
metab.adjP <- predict_metab(adjP.mm)
do.adjP <- predict_DO(adjP.mm)

param.warnings<-get_params(adjP.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(adjP.mm) %>%
  select(date,warnings,errors) 



#Use Values from adjP for new daily estimates
daily.name <- mm_name(type="mle",GPP_fun = "satlight", ER_fun = "q10temp")
daily.specs <- specs(daily.name)
daily.mm <- metab(daily.specs,three.dat,data_daily = select(get_params(adjP.mm),date,init.Pmax = Pmax, init.alpha = alpha))

get_params(daily.mm)
plot_metab_preds(daily.mm)
plot_DO_preds(daily.mm)
metab.adjP <- predict_metab(daily.mm)
do.adjP <- predict_DO(daily.mm)

param.warnings<-get_params(daily.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(daily.mm) %>%
  select(date,warnings,errors) 



#boring mle
boring.name <- mm_name(type="mle")
boring.specs <- specs(boring.name)
boring.mm <- metab(boring.specs,three.dat)

get_params(boring.mm)
plot_metab_preds(boring.mm)
plot_DO_preds(boring.mm)
metab.adjP <- predict_metab(boring.mm)
do.adjP <- predict_DO(boring.mm)

param.warnings<-get_params(boring.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(boring.mm) %>%
  select(date,warnings,errors)


#spec.adj.boring mle
spec.adj.name <- mm_name(type="mle")
spec.adj.specs <- specs(spec.adj.name,init.GPP.daily = 2,init.ER.daily = -1,init.K600.daily = 3)
spec.adj.mm <- metab(spec.adj.specs,three.dat)

get_params(spec.adj.mm)
plot_metab_preds(spec.adj.mm)
plot_DO_preds(spec.adj.mm)
metab.adjP <- predict_metab(spec.adj.mm)
do.adjP <- predict_DO(spec.adj.mm)

param.warnings<-get_params(spec.adj.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(spec.adj.mm) %>%
  select(date,warnings,errors)



#boring q10 model
resp.name <- mm_name(type="mle",ER_fun = "q10temp")
resp.specs <- specs(resp.name)
resp.mm <- metab(resp.specs,three.dat)

get_params(resp.mm)
plot_metab_preds(resp.mm)
plot_DO_preds(resp.mm)
metab.adjP <- predict_metab(resp.mm)
do.adjP <- predict_DO(resp.mm)

param.warnings<-get_params(resp.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(resp.mm) %>%
  select(date,warnings,errors)


#daily boring resp model
daily.resp.name <- mm_name(type="mle",ER_fun = "q10temp")
daily.resp.specs <- specs(daily.resp.name)
daily.resp.mm <- metab(daily.resp.specs,three.dat,data_daily = select(get_params(resp.mm),date,init.ER20 = ER20))

get_params(daily.resp.mm)
plot_metab_preds(daily.resp.mm)
plot_DO_preds(daily.resp.mm)
metab.adjP <- predict_metab(daily.resp.mm)
do.adjP <- predict_DO(daily.resp.mm)

param.warnings<-get_params(daily.resp.mm) %>%
  select(date,warnings,errors) #warnings about Pmax and alpha, satlight params
metab.warnings <- predict_metab(daily.resp.mm) %>%
  select(date,warnings,errors)
