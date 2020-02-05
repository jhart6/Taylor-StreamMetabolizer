#Dissolved Oxygen sensor comparison
#4 minidots: mouth, garden, hound, and culvert from Taylor Creek
#1 EXO2
#Sensors synced at 5min recording frequency & placed in bucket of tap water for >3 days
#5 February 2020
#JAH

library(dplyr)
library(tidyr)
library(streamMetabolizer)

setwd("~/Dropbox/GitHub Repos/Taylor-StreamMetabolizer/DO Bucket Tests/200131_DOBucketTest/")
exo <- read.csv("exo_concat.csv")
culvert <- read.csv("culvert_concat.csv")
garden <- read.csv("garden_concat.csv")
hound <- read.csv("hound_concat.csv")
mouth <- read.csv("mouth_concat.csv")
do.comp <- exo %>%
  select(DateTime,exo2.sat,exo2.mgL)

garden.mgL <- garden$DO..mg.l.
garden.temp <- garden$T..deg.C.
garden.eqsat <- calc_DO_sat(garden.temp,press = 1013,sal = 0)
garden.sat <- (garden.mgL/garden.eqsat)*100

culvert.mgL <- culvert$DO..mg.l.
culvert.temp <- culvert$T..deg.C.
culvert.eqsat <- calc_DO_sat(culvert.temp,press=1013,sal=0)
culvert.sat <- (culvert.mgL/culvert.eqsat) * 100

hound.mgL <- hound$DO..mg.l.
hound.temp <- hound$T..deg.C.
hound.eqsat <- calc_DO_sat(hound.temp,press = 1013, sal=0)
hound.sat <- (hound.mgL/hound.eqsat) *100

mouth.mgL <- mouth$DO..mg.l.
mouth.temp <- mouth$T..deg.C.
mouth.eqsat <- calc_DO_sat(mouth.temp,press=1013, sal=0)
mouth.sat <- (mouth.mgL/mouth.eqsat) * 100

do.comp <- cbind(do.comp,garden.sat,garden.mgL,culvert.sat,culvert.mgL,hound.sat,hound.mgL,mouth.sat,mouth.mgL)


#####Plots#####
#time series plots
quartz()
plot(do.comp$DateTime,do.comp$exo2.sat,ylim=c(60,120),xlab="DateTime",ylab="DO %Sat")
points(do.comp$DateTime,do.comp$garden.sat,col="green")
points(do.comp$DateTime,do.comp$culvert.sat,col="orange")
points(do.comp$DateTime,do.comp$hound.sat,col="blue")
points(do.comp$DateTime,do.comp$mouth.sat,col="red")
legend("bottomright",c("EXO2","Garden MiniDOT","Culvert MiniDOT","Hound MiniDOT","Mouth MiniDOT"),lty=c(1,1,1,1,1),col=c('black',"green","orange","blue","red"))

quartz()
plot(do.comp$DateTime,do.comp$exo2.mgL,ylim=c(2,12),xlab="DateTime",ylab="DO mg/L")
points(do.comp$DateTime,do.comp$garden.mgL,col="green")
points(do.comp$DateTime,do.comp$culvert.mgL,col="orange")
points(do.comp$DateTime,do.comp$hound.mgL,col="blue")
points(do.comp$DateTime,do.comp$mouth.mgL,col="red")
legend("bottomright",c("EXO2","Garden MiniDOT","Culvert MiniDOT","Hound MiniDOT","Mouth MiniDOT"),lty=c(1,1,1,1,1),col=c('black',"green","orange","blue","red"))

#calibration comparison
quartz()
plot(do.comp$exo2.sat,do.comp$garden.sat,ylim=c(75,120),xlim=c(75,120),xlab="EXO2 %Sat",ylab="%Sat",col="green")
abline(0,1,lty=1,cex=2)
points(do.comp$exo2.sat,do.comp$culvert.sat,col="orange")
points(do.comp$exo2.sat,do.comp$hound.sat,col="blue")
points(do.comp$exo2.sat,do.comp$mouth.sat,col="red")
legend("topleft",c("Garden MiniDOT 4.65%","Culvert MiniDOT 5.97%","Hound MiniDOT 8.79%","Mouth MiniDOT 13.85%"),lty=c(1,1,1,1),col=c("green","orange","blue","red"))

quartz()
plot(do.comp$exo2.mgL,do.comp$garden.mgL,ylim=c(7,11),xlim=c(7,11),xlab="EXO2 mg/L",ylab="mg/L",col="green")
abline(0,1,lty=1,cex=2)
points(do.comp$exo2.mgL,do.comp$culvert.mgL,col="orange")
points(do.comp$exo2.mgL,do.comp$hound.mgL,col="blue")
points(do.comp$exo2.mgL,do.comp$mouth.mgL,col="red")
legend("topleft",c("Garden MiniDOT 4.74%","Culvert MiniDOT 6.08%","Hound MiniDOT 8.79%","Mouth MiniDOT 13.75%"),lty=c(1,1,1,1),col=c("green","orange","blue","red"))


####Estimate Offset####
garden.mgL.offset <- mean(((do.comp$exo2.mgL-do.comp$garden.mgL)/do.comp$exo2.mgL) * 100) #4.746782
culvert.mgL.offset <- mean(((do.comp$exo2.mgL-do.comp$culvert.mgL)/do.comp$exo2.mgL) * 100) #6.084556
hound.mgL.offset <- mean(((do.comp$exo2.mgL-do.comp$hound.mgL)/do.comp$exo2.mgL) * 100) #8.796492
mouth.mgL.offset <- mean(((do.comp$exo2.mgL-do.comp$mouth.mgL)/do.comp$exo2.mgL) * 100) #13.85024

garden.sat.offset <- mean(((do.comp$exo2.sat-do.comp$garden.sat)/do.comp$exo2.sat) * 100) #4.658371
culvert.sat.offset <- mean(((do.comp$exo2.sat-do.comp$culvert.sat)/do.comp$exo2.sat) * 100) #5.979058
hound.sat.offset <- mean(((do.comp$exo2.sat-do.comp$hound.sat)/do.comp$exo2.sat) * 100) #8.794387
mouth.sat.offset <- mean(((do.comp$exo2.sat-do.comp$mouth.sat)/do.comp$exo2.sat) * 100) #13.75172

####Adjust Observed Data####
garden.mgL.adj <- do.comp$garden.mgL * (1+(garden.mgL.offset/100))
culvert.mgL.adj <- do.comp$culvert.mgL * (1 +(culvert.mgL.offset/100))
hound.mgL.adj <- do.comp$hound.mgL * (1+(hound.mgL.offset/100))
mouth.mgL.adj <- do.comp$mouth.mgL * (1+(mouth.mgL.offset/100))

garden.sat.adj <- do.comp$garden.sat * (1+(garden.sat.offset/100))
culvert.sat.adj <- do.comp$culvert.sat * (1+(culvert.sat.offset/100))
hound.sat.adj <- do.comp$hound.sat * (1+(hound.sat.offset/100))
mouth.sat.adj <- do.comp$mouth.sat * (1+(mouth.sat.offset/100))

#####Plot Adjusted Data####
quartz()
plot(do.comp$DateTime,do.comp$exo2.sat,ylim=c(60,120),xlab="DateTime",ylab="DO %Sat")
points(do.comp$DateTime,garden.sat.adj,col="green")
points(do.comp$DateTime,culvert.sat.adj,col="orange")
points(do.comp$DateTime,hound.sat.adj,col="blue")
points(do.comp$DateTime,mouth.sat.adj,col="red")
legend("bottomright",c("EXO2","Garden MiniDOT Adj","Culvert MiniDOT Adj","Hound MiniDOT Adj","Mouth MiniDOT Adj"),lty=c(1,1,1,1,1),col=c('black',"green","orange","blue","red"))

quartz()
plot(do.comp$DateTime,do.comp$exo2.mgL,ylim=c(2,12),xlab="DateTime",ylab="DO mg/L")
points(do.comp$DateTime,garden.mgL.adj,col="green")
points(do.comp$DateTime,culvert.mgL.adj,col="orange")
points(do.comp$DateTime,hound.mgL.adj,col="blue")
points(do.comp$DateTime,mouth.mgL.adj,col="red")
legend("bottomright",c("EXO2","Garden MiniDOT Adj","Culvert MiniDOT Adj","Hound MiniDOT Adj","Mouth MiniDOT Adj"),lty=c(1,1,1,1,1),col=c('black',"green","orange","blue","red"))

#calibration comparison
quartz()
plot(do.comp$exo2.sat,garden.sat.adj,ylim=c(85,120),xlim=c(85,120),xlab="EXO2 %Sat",ylab="%Sat",col="green")
abline(0,1,lty=1,cex=2)
points(do.comp$exo2.sat,culvert.sat.adj,col="orange")
points(do.comp$exo2.sat,hound.sat.adj,col="blue")
points(do.comp$exo2.sat,mouth.sat.adj,col="red")
legend("topleft",c("Garden MiniDOT Adj","Culvert MiniDOT Adj","Hound MiniDOT Adj","Mouth MiniDOT Adj"),lty=c(1,1,1,1),col=c("green","orange","blue","red"))

quartz()
plot(do.comp$exo2.mgL,garden.mgL.adj,ylim=c(7,11),xlim=c(7,11),xlab="EXO2 mg/L",ylab="mg/L",col="green")
abline(0,1,lty=1,cex=2)
points(do.comp$exo2.mgL,culvert.mgL.adj,col="orange")
points(do.comp$exo2.mgL,hound.mgL.adj,col="blue")
points(do.comp$exo2.mgL,mouth.mgL.adj,col="red")
legend("topleft",c("Garden MiniDOT Adj","Culvert MiniDOT Adj","Hound MiniDOT Adj","Mouth MiniDOT Adj"),lty=c(1,1,1,1),col=c("green","orange","blue","red"))

