library(Rmisc)
library(dplyr)
library(tidyr)
library(lubridate)
library(forecast)
library(tseries)
library(plyr)
library(VGAM)
library(pscl)
library(ggpointdensity)
library(MASS)
library(gridExtra)
library(mgcv)
library(tidyverse)
library(forecast)
library(tsibble)
library(tscount)



set.seed(3)
bikedata<-read.csv("../Real-data/hour.csv")


bikedata<-bikedata %>% 
  filter(yr==1)
bikedata<-bikedata[,!grepl("instant",colnames(bikedata))]
bikedata<-bikedata[,!grepl("casual",colnames(bikedata))]
bikedata<-bikedata[,!grepl("registered",colnames(bikedata))]
bikedata<-bikedata[,!grepl("atemp",colnames(bikedata))]
bikedata<-bikedata[,!grepl("yr",colnames(bikedata))]
bikedata<-bikedata[,!grepl("holiday",colnames(bikedata))]


bikedata$dteday <- as.Date(bikedata$dteday, format = "%Y-%m-%d")
bikedata$winter<-ifelse(bikedata$season==1,1,0)
bikedata<-bikedata[,!grepl("season",colnames(bikedata))]


bikedata$dteday<-as.Date(bikedata$dteday)
bikedata$hr<-as.numeric(bikedata$hr)
bikedata$hour<-bikedata$hr
bikedata <- bikedata %>% mutate(datetime=
                                  as.POSIXct(as.character(paste(bikedata$dteday, bikedata$hr)), 
                                                    format="%Y-%m-%d %H"))
bikedata<-bikedata[,!grepl("hr",colnames(bikedata))]
bikedata<-bikedata[,!grepl("dteday",colnames(bikedata))]
bikedata<-bikedata[,!grepl("mnth",colnames(bikedata))]

xpart<-bikedata %>% select(weekday,workingday,weathersit,temp,
                           windspeed,hum,winter)


##################################################################################
########################## N B Distribution ######################################
##################################################################################

modeltsnb1<-tsglm(bikedata$cnt,distr="nbinom",model=list(past_obs=1),
               xreg = xpart,link = "log")
modeltsnb1$coefficients

i<-10
lineari<-1.269587845+log(bikedata$cnt[i-1]+1)*0.771530373+
  bikedata$weekday[i]*0.003511883+
  bikedata$weathersit[i]*(-0.029266388)+
  bikedata$temp[i]*0.219233355+
  bikedata$windspeed[i]*0.084295495+
  bikedata$hum[i]*(-0.130842873)+
  bikedata$winter[i]*(-0.073513387)
exp(lineari)
bikedata$nbfitted1[i]

modeltsnb2<-tsglm(bikedata$cnt,distr="nbinom",model=list(past_obs=2),
                  xreg = xpart,link = "log")
modeltsnb3<-tsglm(bikedata$cnt,distr="nbinom",model=list(past_obs=3),
                  xreg = xpart,link = "log")
modeltsnb33<-tsglm(bikedata$cnt,distr="nbinom",
                   model=list(past_obs=3,past_mean=3),
                  xreg = xpart)



bikedata$nbfitted1<-modeltsnb1$fitted.values
bikedata$nbfitted2<-modeltsnb2$fitted.values
bikedata$nbfitted3<-modeltsnb3$fitted.values
bikedata$nbfitted33<-modeltsnb33$fitted.values

disp1<-1/modeltsnb1$sigmasq
disp2<-1/modeltsnb2$sigmasq
disp3<-1/modeltsnb3$sigmasq
disp33<-1/modeltsnb33$sigmasq

range_tsnb1<-matrix(NA,nrow = length(bikedata$cnt),ncol = 2)
range_tsnb2<-matrix(NA,nrow = length(bikedata$cnt),ncol = 2)
range_tsnb3<-matrix(NA,nrow = length(bikedata$cnt),ncol = 2)
range_tsnb33<-matrix(NA,nrow = length(bikedata$cnt),ncol = 2)


# ingarch 1,0
for (i in 1:length(bikedata$cnt)) {
  range_tsnb1[i,]<-c(pnbinom(bikedata$cnt[i]-1,
                                 size=disp1,
                                 mu=modeltsnb1$fitted.values[i]),
                         pnbinom(bikedata$cnt[i],
                                 size=disp1,
                                 mu=modeltsnb1$fitted.values[i]))
}
# ingarch 2,0
for (i in 1:length(bikedata$cnt)) {
  range_tsnb2[i,]<-c(pnbinom(bikedata$cnt[i]-1,
                             size=disp2,
                             mu=modeltsnb2$fitted.values[i]),
                     pnbinom(bikedata$cnt[i],
                             size=disp2,
                             mu=modeltsnb2$fitted.values[i]))
}

# ingarch 3,0
for (i in 1:length(bikedata$cnt)) {
  range_tsnb3[i,]<-c(pnbinom(bikedata$cnt[i]-1,
                             size=disp3,
                             mu=modeltsnb3$fitted.values[i]),
                     pnbinom(bikedata$cnt[i],
                             size=disp3,
                             mu=modeltsnb3$fitted.values[i]))
}
# ingarch 3,3
for (i in 1:length(bikedata$cnt)) {
  range_tsnb33[i,]<-c(pnbinom(bikedata$cnt[i]-1,
                             size=disp33,
                             mu=modeltsnb33$fitted.values[i]),
                     pnbinom(bikedata$cnt[i],
                             size=disp33,
                             mu=modeltsnb33$fitted.values[i]))
}



tsnbnumbers1<-matrix(NA,nrow=length(bikedata$cnt),ncol = 11)
for (h in 1:nrow(tsnbnumbers1)) {
  for (a in 1:ncol(tsnbnumbers1)) {
    tsnbnumbers1[h,a]<-range_tsnb1[h,1]+
      (range_tsnb1[h,2]-range_tsnb1[h,1])/10*(a-1)
  }
}

tsnbnumbers2<-matrix(NA,nrow=length(bikedata$cnt),ncol = 11)
for (h in 1:nrow(tsnbnumbers2)) {
  for (a in 1:ncol(tsnbnumbers2)) {
    tsnbnumbers2[h,a]<-range_tsnb2[h,1]+
      (range_tsnb2[h,2]-range_tsnb2[h,1])/10*(a-1)
  }
}

tsnbnumbers3<-matrix(NA,nrow=length(bikedata$cnt),ncol = 11)
for (h in 1:nrow(tsnbnumbers3)) {
  for (a in 1:ncol(tsnbnumbers3)) {
    tsnbnumbers3[h,a]<-range_tsnb3[h,1]+
      (range_tsnb3[h,2]-range_tsnb3[h,1])/10*(a-1)
  }
}

tsnbnumbers33<-matrix(NA,nrow=length(bikedata$cnt),ncol = 11)
for (h in 1:nrow(tsnbnumbers33)) {
  for (a in 1:ncol(tsnbnumbers33)) {
    tsnbnumbers33[h,a]<-range_tsnb33[h,1]+
      (range_tsnb33[h,2]-range_tsnb33[h,1])/10*(a-1)
  }
}

tsnb_vector1<-qnorm(as.vector(tsnbnumbers1))
tsnb_vector2<-qnorm(as.vector(tsnbnumbers2))
tsnb_vector3<-qnorm(as.vector(tsnbnumbers3))
tsnb_vector33<-qnorm(as.vector(tsnbnumbers33))

# for plot purpose
rephour<-rep(bikedata$hour,11)
tsnb_heatmapdata<-cbind.data.frame(rephour,
                                   tsnb_vector1,tsnb_vector2,
                                   tsnb_vector3,tsnb_vector33)


heatmaptsnb1_hour<-ggplot(tsnb_heatmapdata, aes(rephour,tsnb_vector1)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("hour")+ylab("")+
  labs(title = "(a) INGARCH(1,0)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmaptsnb2_hour<-ggplot(tsnb_heatmapdata, aes(rephour,tsnb_vector2)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("hour")+ylab("")+
  labs(title = "(b) INGARCH(2,0)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmaptsnb3_hour<-ggplot(tsnb_heatmapdata, aes(rephour,tsnb_vector3)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("hour")+ylab("")+
  labs(title = "(c) INGARCH(3,0)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmaptsnb33_hour<-ggplot(tsnb_heatmapdata, aes(rephour,tsnb_vector33)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("hour")+ylab("")+
  labs(title = "(d) INGARCH(3,3)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

########################Figure S20#####################################

grid.arrange(heatmaptsnb1_hour,heatmaptsnb2_hour,
             heatmaptsnb3_hour,heatmaptsnb33_hour,ncol=2)









