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

#When a much smaller subset (n = 250) is randomly selected from each of the original data sets. 
#The same analysis is repeated.


set.seed(6)


whitewine<-read.csv("./Real-data/winequality-white.csv",sep = ";")
subsample<-sample(c(1:nrow(whitewine)),size = 250,replace = F)
whitewine<-whitewine[subsample,]


#alcohol+volatile.acidity+residual.sugar+free.sulfur.dioxide+density+pH+sulphates+fixed.acidity+citric.acid
model1ww<- vglm(quality~volatile.acidity+
                alcohol+sulphates+fixed.acidity+
                residual.sugar+free.sulfur.dioxide+
                pH+density,
              family=acat(reverse=TRUE, parallel=TRUE),data =whitewine)

fr1_sub_ww<-fresiduals(model1ww)
#######Function residuals-vs-covariate plot


heatmap6_norm<-fresplot(fr1_sub_ww,whitewine$free.sulfur.dioxide,
                        scale ="normal",
                        xl=0,xp=150,heatmapcut = 11,yl=-3,yp=3,
                        title = "(a) free.sulfur.dioxide",xlabs = "")
  

bikedata<-read.csv("./Real-data/hour.csv")


bikedata<-bikedata %>% 
  filter(yr==1)

bikedata$winter<-ifelse(bikedata$season==1,1,0)

subsample<-sample(c(1:nrow(bikedata)),size = 250,replace = F)

bikedata<-bikedata[subsample,]

model1bike<-glm(cnt~hr+workingday+weathersit+temp+hum+windspeed+winter,family = "poisson",data = bikedata)

### Functional residuals

fr_sub_bike<-fresiduals(model1bike)


######heatmap


heatmap2_norm<-fresplot(fr_sub_bike,bikedata$hr,
                        scale ="normal",
                        xl=0,xp=24,heatmapcut = 11,yl=-20,yp=20,
                        title = "(b) hour",xlabs = "")


##########additive model


model_gam<-gam(cnt~winter+s(hr)+workingday+weathersit+
                 s(temp)+s(hum)+s(windspeed),
               family = poisson,
               data = bikedata)


summary(model_gam)

### Functional residuals

fr_sub_gam_bike<-fresiduals(model_gam)


###########

ff<-ffplot(fr_sub_gam_bike,title = "(c) Fn-Fn plot")

#####quasipoisson


model_gam_quasi<-gam(cnt~winter+s(hr)+workingday+weathersit+
                       s(temp)+s(hum)+s(windspeed),
                     family = quasipoisson,
                     data = bikedata)
summary(model_gam_quasi)

fr_sub_gam_quasi<-fresiduals(model_gam_quasi)

heatmap2_norm_gam_quasi<-fresplot(fr_sub_gam_quasi,bikedata$hr,
                                  scale ="normal",
                                  xl=0,xp=24,heatmapcut = 11,yl=-3,yp=3,
                                  title = "(d) hour",xlabs = "")

################################################
##Figure S24 ###################################
################################################
grid.arrange(heatmap6_norm,heatmap2_norm,
             ff,heatmap2_norm_gam_quasi,nrow=2)
