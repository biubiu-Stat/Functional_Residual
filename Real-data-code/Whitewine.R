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
##########################################################################
##########################################################################
######################White Wine Data Analysis############################
##########################################################################
##########################################################################
set.seed(3)

whitewine<-read.csv("./Real-data/winequality-white.csv",sep = ";")

model1<- vglm(quality~volatile.acidity+
                alcohol+sulphates+fixed.acidity+
                residual.sugar+free.sulfur.dioxide+
                pH+density,
              family=acat(reverse=TRUE, parallel=TRUE),data =whitewine)

################################Initial Model for Table S1######################
summary(model1)

################################Functional Residual for Initial Model###########
fr1_whitewine<-fresiduals(model1)

#######heatmap part

heatmap2_norm<-fresplot(fr1_whitewine,whitewine$fixed.acidity,
                        scale ="normal",
                        xl=4,xp=15,heatmapcut = 11,
                        title = "(a) fixed.acidity",xlabs = "")


heatmap5_norm<-fresplot(fr1_whitewine,whitewine$residual.sugar,
                        scale ="normal",
                        xl=0,xp=70,heatmapcut = 11,yl=-4,yp=4,
                        title = "(b) residual.sugar",xlabs = "")

min(whitewine$density)

max(whitewine$density)
heatmap8_norm<-fresplot(fr1_whitewine,whitewine$density,
                        scale ="normal",
                        xl=0.985,xp=1.04,heatmapcut = 11,yl=-4,yp=4,
                        title = "(c) density",xlabs = "")

##########################################################################
###############################Figure S11#################################
##########################################################################

par(mfrow = c(2, 2))

boxplot(whitewine$fixed.acidity, main="(a) fixed.acidity" , horizontal = TRUE,col="#69b3a2", boxwex=0.4)
stripchart(c(max(whitewine$fixed.acidity),11.8),pch = 4, col ="red", vertical = FALSE, add = TRUE,cex=5)

boxplot(whitewine$residual.sugar,horizontal = TRUE, main="(b) residual.sugar" , col="#69b3a2", boxwex=0.4)
stripchart(max(whitewine$residual.sugar),pch = 4, col ="red", vertical = FALSE, add = TRUE,cex=5)
boxplot(whitewine$density, main="(c) density" ,horizontal = TRUE, col="#69b3a2", boxwex=0.4 , main="")
stripchart(c(1.0103,1.03898),pch = 4, col ="red", vertical = FALSE, add = TRUE,cex=5)

boxplot(whitewine$free.sulfur.dioxide, main="(d) free.sulfur.dioxide" ,horizontal = TRUE, col="#69b3a2", boxwex=0.4 , main="")

# Restore default layout
par(mfrow = c(1, 1))




#############delete outlier

whitewine_rmout<-cbind.data.frame(fr1_whitewine,whitewine) %>%
  filter(residual.sugar<=50) %>%
  filter(density<1.01) %>%
  filter(fixed.acidity<11)
fr1_update<-as.matrix(whitewine_rmout[,1:2])

heatmap2_norm_v2<-fresplot(fr1_update,whitewine_rmout$fixed.acidity,
                           scale ="normal",
                           xl=4,xp=11,heatmapcut = 11,
                           title = "(a*) fixed.acidity",xlabs = "")
  
heatmap5_norm_v2<-fresplot(fr1_update,whitewine_rmout$residual.sugar,
                           scale ="normal",
                           xl=0,xp=24,heatmapcut = 11,yl=-4,yp=4,
                           title = "(b*) residual.sugar",xlabs = "")

heatmap8_norm_v2<-fresplot(fr1_update,whitewine_rmout$density,
                           scale ="normal",yl=-4,yp=4,
                           xl=0.986,xp=1.003,heatmapcut = 11,
                           title = "(c*) density",xlabs = "")

####################################################
#####################Figure S10#####################
####################################################

grid.arrange(heatmap2_norm,heatmap5_norm,heatmap8_norm,
             heatmap2_norm_v2,heatmap5_norm_v2,heatmap8_norm_v2,nrow=2)



#################### Add square term #######################

whitewine_rmout$free.sulfur.dioxide2<-whitewine_rmout$free.sulfur.dioxide^2

model2<- vglm(quality~volatile.acidity+
                alcohol+sulphates+fixed.acidity+
                residual.sugar+free.sulfur.dioxide+
                pH+density+free.sulfur.dioxide2,
              family=acat(reverse=TRUE, parallel=TRUE),data =whitewine_rmout)

################################Final Model for Table S1######################

summary(model2)

##try one by one AIC=11023.12
#chlorides+volatile.acidity+total.sulfur.dioxide+alcohol
#sulphates+fixed.acidity+citric.acid+residual.sugar+free.sulfur.dioxide+pH+density



################################Functional Residual for Updated Model###########
fr2_whitewine<-fresiduals(model2)

heatmapbefore<-fresplot(fr1_update,whitewine_rmout$free.sulfur.dioxide,
                        scale ="normal",
                        xl=0,xp=300,heatmapcut = 11,yl=-5,yp=5,
                        title = "(a) before",xlabs = "")



heatmap_norm_after_quard<-fresplot(fr2_whitewine,whitewine_rmout$free.sulfur.dioxide,
                                   scale ="normal",
                                   xl=0,xp=300,heatmapcut = 11,yl=-5,yp=5,
                                   title = "(b) after",xlabs = "")



####################################################
#####################Figure S12#####################
####################################################
multiplot(heatmapbefore,heatmap_norm_after_quard,cols = 2) 

