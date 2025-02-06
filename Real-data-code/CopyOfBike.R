library(Rmisc)
library(plyr)
library(ggplot2)
library(VGAM)
library(pscl)
library(ggpointdensity)
library(MASS)
library(dplyr)
library(gridExtra)
library(mgcv)

set.seed(3)
bikedata<-read.csv(here::here("./Real-data/hour.csv"))

### Select the data in 2012
bikedata<-bikedata %>% 
  filter(yr==1)


# Remove columns related to casual users and registered users, as they are not needed for this analysis
bikedata <- bikedata[, !grepl("casual", colnames(bikedata))]
bikedata <- bikedata[, !grepl("registered", colnames(bikedata))]

# Remove the "atemp" column to avoid redundancy of temp
bikedata <- bikedata[, !grepl("atemp", colnames(bikedata))]

# Remove the "yr" column (year), we select the data in 2012
bikedata <- bikedata[, !grepl("yr", colnames(bikedata))]

# Remove the "dteday" column (date), as they are not needed for this analysis
bikedata <- bikedata[, !grepl("dteday", colnames(bikedata))]

# Remove the "holiday" column, as holidays may not be considered for this analysis
bikedata <- bikedata[, !grepl("holiday", colnames(bikedata))]

# Remove the "instant" column
bikedata <- bikedata[, !grepl("instant", colnames(bikedata))]

# Remove the "mnth" column (month)
bikedata <- bikedata[, !grepl("mnth", colnames(bikedata))]

# Remove the "weekday" column, we have workday
bikedata <- bikedata[, !grepl("weekday", colnames(bikedata))]
### Create the winter variable
bikedata$winter<-ifelse(bikedata$season==1,1,0)
bikedata<-bikedata[,!grepl("season",colnames(bikedata))]

### Initial model with all the variables
model1<-glm(cnt~.,family = "poisson",data = bikedata)
fr1 <- fresiduals(model1)
#####Initial model for Table S3 #####
summary(model1)
fittedy1<-model1$fitted.values #predicted value
n<-length(fittedy1)

#### Range for functional residuals
bikerange<-matrix(NA,nrow=n,ncol = 2)
for (i in 1:n) {
  bikerange[i,]<-c(ppois(bikedata$cnt[i]-1,fittedy1[i]),
                   ppois(bikedata$cnt[i],fittedy1[i]))
}

#### cut each range for plot
pre<-matrix(NA,nrow=n,ncol = 11)
for (h in 1:n) {
  for (a in 1:ncol(pre)) {
    pre[h,a]<-bikerange[h,1]+(bikerange[h,2]-bikerange[h,1])/10*(a-1)
  }
}
qnumbervector<-qnorm(as.vector(pre))
### add random noise for visualization 
repwinter<-rep(bikedata$winter,11)+runif(length(qnumbervector), 0, 0.01)
### repeat variables for visualization 
rephr<-rep(bikedata$hr,11)
repworkingday<-rep(bikedata$workingday,11)
repweathersit<-rep(bikedata$weathersit,11)
reptemp<-rep(bikedata$temp,11)
rephum<-rep(bikedata$hum,11)
repwindspeed<-rep(bikedata$windspeed,11)


bikeheatmapdata<-cbind.data.frame(repwinter,
                                  repworkingday,repweathersit,reptemp,rephr,
                                  rephum,repwindspeed,
                                  qnumbervector)

#bikeheatmapdata$repwinter<-bikeheatmapdata$repwinter+runif(nrow(bikeheatmapdata), -0.05, 0.05)
heatmap_winter<-fresplot(fr1,bikedata$winter,
                         title = "(e)winter",scale = "normal",xl=0,xp=1,
                         xlabs = "", heatmapcut=11)

heatmap_winter<-ggplot(bikeheatmapdata, aes(repwinter,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(e)winter")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
heatmap_hour<-ggplot(bikeheatmapdata, aes(rephr,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(a) hour")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_workingday<-ggplot(bikeheatmapdata, aes(repworkingday,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(f) workingday")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap4_norm<-ggplot(bikeheatmapdata, aes(repweathersit,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("weathersit")+ylab("")+
  labs(title = "Functional residuals")

heatmap_temp<-ggplot(bikeheatmapdata, aes(reptemp,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(b) temp")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_humidity<-ggplot(bikeheatmapdata, aes(rephum,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(c) humidity")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_windspeed<-ggplot(bikeheatmapdata, aes(repwindspeed,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(d) windspeed")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

###############################################################################
###############################Figure S13 #####################################
###############################################################################
grid.arrange(heatmap_hour,heatmap_temp,heatmap_humidity,heatmap_windspeed,
             heatmap_winter,
             heatmap_workingday,ncol=2)

#####################Fn-Fn plot for initial model##############

t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-as.data.frame(bikerange)
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest$V1[h],max=range1fortest$V2[h])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1_initial<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(a) Initial model")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


########## Additive Poisson model


model_gam<-gam(cnt~winter+s(hr)+workingday+weathersit+
              s(temp)+s(hum)+s(windspeed),
            family = poisson,
            data = bikedata)

summary(model_gam)
fittedy_gam<-model_gam$fitted.values

n<-length(fittedy_gam)



bikerange_gam<-matrix(NA,nrow=n,ncol = 2)
for (i in 1:n) {
  bikerange_gam[i,]<-c(ppois(bikedata$cnt[i]-1,fittedy_gam[i]),
                   ppois(bikedata$cnt[i],fittedy_gam[i]))
}

pre_gam<-matrix(NA,nrow=n,ncol = 11)
for (h in 1:n) {
  for (a in 1:ncol(pre_gam)) {
    pre_gam[h,a]<-bikerange_gam[h,1]+(bikerange_gam[h,2]-bikerange_gam[h,1])/10*(a-1)
  }
}
qnumbervector_gam<-qnorm(as.vector(pre_gam))

bikeheatmapdata_gam<-cbind.data.frame(repwinter,
                                  repworkingday,repweathersit,reptemp,rephr,
                                  rephum,repwindspeed,
                                  qnumbervector_gam)
heatmap1_norm_gam<-ggplot(bikeheatmapdata_gam, aes(repwinter,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("Winter")+ylab("")+
  labs(title = "Functional residuals")


heatmap2_norm_gam<-ggplot(bikeheatmapdata_gam, aes(rephr,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(a) hour")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap3_norm_gam<-ggplot(bikeheatmapdata_gam, aes(repworkingday,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("workingday")+ylab("")+
  labs(title = "Functional residuals")

heatmap4_norm_gam<-ggplot(bikeheatmapdata_gam, aes(repweathersit,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("weathersit")+ylab("")+
  labs(title = "Functional residuals")

heatmap5_norm_gam<-ggplot(bikeheatmapdata_gam, aes(reptemp,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(b) temp")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
heatmap6_norm_gam<-ggplot(bikeheatmapdata_gam, aes(rephum,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(c) humidity")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap7_norm_gam<-ggplot(bikeheatmapdata_gam, aes(repwindspeed,qnumbervector_gam)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(d) windspeed")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

###############################################################################
###############################Figure S15 #####################################
###############################################################################

grid.arrange(heatmap2_norm_gam,heatmap5_norm_gam,heatmap6_norm_gam,
             heatmap7_norm_gam,nrow=2)
########### Fn-Fn for Intermediate#########

t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-as.data.frame(bikerange_gam)
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest$V1[h],max=range1fortest$V2[h])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1_gam<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(b) Intermediate model")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))




#####quasipoisson


model_gam_quasi<-gam(cnt~winter+s(hr)+workingday+weathersit+
                 s(temp)+s(hum)+s(windspeed),
               family = quasipoisson,
               data = bikedata)

############### Final model for Table S3 ###############
summary(model_gam_quasi)

gamfitted<-model_gam_quasi$fitted.values

s_gam_quasi<-summary(model_gam_quasi)
range_gam_quasi<-matrix(NA,nrow = length(model_gam_quasi$fitted.values),ncol = 2)
for (i in 1:length(model_gam_quasi$fitted.values)) {
  range_gam_quasi[i,]<-c(pnbinom(bikedata$cnt[i]-1,
                                 size=model_gam_quasi$fitted.values[i]/s_gam_quasi$dispersion,
                                 mu=model_gam_quasi$fitted.values[i]),
                         pnbinom(bikedata$cnt[i],
                                 size=model_gam_quasi$fitted.values[i]/s_gam_quasi$dispersion,
                                 mu=model_gam_quasi$fitted.values[i]))
}

gam_quasi_numbers1<-matrix(NA,nrow=length(model_gam_quasi$fitted.values),ncol = 11)
for (h in 1:nrow(gam_quasi_numbers1)) {
  for (a in 1:ncol(gam_quasi_numbers1)) {
    gam_quasi_numbers1[h,a]<-range_gam_quasi[h,1]+(range_gam_quasi[h,2]-range_gam_quasi[h,1])/10*(a-1)
  }
}

gam_quasi_vector<-qnorm(as.vector(gam_quasi_numbers1))
gam_quasi_heatmapdata<-cbind.data.frame(repwinter,
                                  repworkingday,repweathersit,reptemp,rephr,
                                  rephum,repwindspeed,
                                  gam_quasi_vector)

heatmap_winter_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(repwinter,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+xlim(0,1.01)+
  labs(title = "(e) winter")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_hour_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(rephr,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(a) hour")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


heatmap_workingday_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(repworkingday,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(f) workingday")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(repweathersit,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("weathersit")+ylab("")+
  labs(title = "Functional residuals")

heatmap_temp_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(reptemp,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(b) temp")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
heatmap_humidity_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(rephum,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(c) humidity")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_windspeed_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(repwindspeed,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(d) windspeed")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
########################################################################
####################################Figure S16 #########################
########################################################################

grid.arrange(heatmap_hour_gam_quasi,heatmap_temp_gam_quasi,
             heatmap_humidity_gam_quasi,heatmap_windspeed_gam_quasi,
             heatmap_winter_gam_quasi,
             heatmap_workingday_gam_quasi,ncol=2)


######Fn-Fn



t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-as.data.frame(range_gam_quasi)
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest$V1[h],max=range1fortest$V2[h])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1_final<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(c) Final model")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

####################################Figure S14 #################################
grid.arrange(tpoints1_initial,tpoints1_gam,tpoints1_final,nrow=1)


