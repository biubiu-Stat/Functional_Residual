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
source("./functions/ffplot.R")
source("./functions/fresiduals.R")
source("./functions/fresplot.R")


set.seed(3)
bikedata<-read.csv(here::here("./Real-data/hour.csv"))

### Our analysis focuses on the 2012 data (Selecting the 2012 data).
bikedata<-bikedata %>% 
  filter(yr==1)

### Create the winter variable
bikedata$winter<-ifelse(bikedata$season==1,1,0)


#1. information is contained/absorbed

### Initial model with all the variables
model1<-glm(cnt~winter+hr+workingday+weathersit+
              temp+hum+windspeed,family = "poisson",data = bikedata) # the remaining variables are....
## Functional residuals
fr1 <- fresiduals(model1)
#####Initial model for Table S3 #####
summary(model1)

heatmap_winter<-fresplot(fr1,bikedata$winter,
                         title = "(e) winter",scale = "normal",xl=0,xp=1,
                         xlabs = "", heatmapcut=11)

heatmap_hour<-fresplot(fr1, bikedata$hr,
                       title = "(a) hour",scale = "normal",xl=0, xp=24,
                       xlabs = "", heatmapcut=11)

heatmap_workingday<-fresplot(fr1, bikedata$workingday,
                             title = "(f) workingday",scale = "normal",
                             xl=0, xp=1,
                             xlabs = "", heatmapcut=11)

heatmap_temp<-fresplot(fr1, bikedata$temp,
                       title = "(b) temp",scale = "normal",
                       xl=0, xp=1,
                       xlabs = "", heatmapcut=11)

heatmap_humidity<-fresplot(fr1, bikedata$hum,
                           title = "(c) humidity",scale = "normal",
                           xl=0.16, xp=1,
                           xlabs = "", heatmapcut=11)

heatmap_windspeed<-fresplot(fr1, bikedata$windspeed,
                            title = "(d) windspeed",scale = "normal",
                            xl=0, xp=0.8,
                            xlabs = "", heatmapcut=11)

###############################################################################
############### Figure S13 Functional-residual-vs-covariate plots 
############### for the initial Poisson model fitted to the Captial Bikeshare dataset.
###############################################################################
grid.arrange(heatmap_hour,heatmap_temp,heatmap_humidity,heatmap_windspeed,
             heatmap_winter,
             heatmap_workingday,ncol=2)

#####################Fn-Fn plot for initial model##############

ff1<-ffplot(fr1,title = "(a) Initial model")


########## Additive Poisson model


model_gam<-gam(cnt~winter+s(hr)+workingday+weathersit+
              s(temp)+s(hum)+s(windspeed),
            family = poisson,
            data = bikedata)

######### Functional residuals for Additive Poisson model

fr2<-fresiduals(model_gam)




heatmap2_norm_gam<-fresplot(fr2, bikedata$hr,
                            title = "(a) hour",scale = "normal",xl=0, xp=24,
                            xlabs = "", heatmapcut=11)



heatmap5_norm_gam<-fresplot(fr2, bikedata$temp,
                            title = "(b) temp",scale = "normal",
                            xl=0, xp=1,
                            xlabs = "", heatmapcut=11)

heatmap6_norm_gam<-fresplot(fr2, bikedata$hum,
                            title = "(c) humidity",scale = "normal",
                            xl=0.16, xp=1,
                            xlabs = "", heatmapcut=11)

heatmap7_norm_gam<-fresplot(fr2, bikedata$windspeed,
                            title = "(d) windspeed",scale = "normal",
                            xl=0, xp=0.8,
                            xlabs = "", heatmapcut=11)
###############################################################################
###############################Figure S15 Functional-residual-vs-covariate plots after adding the smoothing functions#####################################
###############################################################################

grid.arrange(heatmap2_norm_gam,heatmap5_norm_gam,heatmap6_norm_gam,
             heatmap7_norm_gam,nrow=2)
########### Fn-Fn for Intermediate#########
ff2<-ffplot(fr2,title = "(b) Intermediate model")

#####quasipoisson


model_gam_quasi<-gam(cnt~winter+s(hr)+workingday+weathersit+
                 s(temp)+s(hum)+s(windspeed),
               family = quasipoisson,
               data = bikedata)
fr3<- fresiduals(model_gam_quasi)

############### Final model for Table S3 ###############
summary(model_gam_quasi)
set.seed(3)
heatmap_winter_gam_quasi<-fresplot(fr3,bikedata$winter ,
                         title = "(e) winter",scale = "normal",xl=0,xp=1.01,
                         xlabs = "", heatmapcut=11,is.binary=TRUE)


heatmap_hour_gam_quasi<-fresplot(fr3, bikedata$hr,
                       title = "(a) hour",scale = "normal",xl=0, xp=24,
                       xlabs = "", heatmapcut=11)


heatmap_workingday_gam_quasi<-fresplot(fr3, bikedata$workingday,
                             title = "(f) workingday",scale = "normal",
                             xl=0, xp=1,
                             xlabs = "", heatmapcut=11)


heatmap_temp_gam_quasi<-fresplot(fr3, bikedata$temp,
                       title = "(b) temp",scale = "normal",
                       xl=0, xp=1,
                       xlabs = "", heatmapcut=11)

heatmap_humidity_gam_quasi<-fresplot(fr3, bikedata$hum,
                           title = "(c) humidity",scale = "normal",
                           xl=0.16, xp=1,
                           xlabs = "", heatmapcut=11)


heatmap_windspeed_gam_quasi<-fresplot(fr3, bikedata$windspeed,
                            title = "(d) windspeed",scale = "normal",
                            xl=0, xp=0.8,
                            xlabs = "", heatmapcut=11)




########################################################################
####################################Figure S16 Functional-residual-vs-covariate plots for the final model.#########################
########################################################################

grid.arrange(heatmap_hour_gam_quasi,heatmap_temp_gam_quasi,
             heatmap_humidity_gam_quasi,heatmap_windspeed_gam_quasi,
             heatmap_winter_gam_quasi,
             heatmap_workingday_gam_quasi,ncol=2)


######Fn-Fn

ff3<-ffplot(fr3,title="(c) Final model")


####################################Figure S14 The Fn-Fn plots for the initial, intermediate, and final models#################################
grid.arrange(ff1,ff2,ff3,nrow=1)


