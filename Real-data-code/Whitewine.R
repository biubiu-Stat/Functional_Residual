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
##########################################################################
##########################################################################
###################################build model############################
##########################################################################
##########################################################################
set.seed(3)

whitewine<-read.csv("../Real-data/winequality-white.csv",sep = ";")

#alcohol+volatile.acidity+residual.sugar+free.sulfur.dioxide+density+pH+sulphates+fixed.acidity+citric.acid
model1<- vglm(quality~volatile.acidity+
                alcohol+sulphates+fixed.acidity+
                residual.sugar+free.sulfur.dioxide+
                pH+density,
              family=acat(reverse=TRUE, parallel=TRUE),data =whitewine)

#alcohol+volatile.acidity+residual.sugar+free.sulfur.dioxide+density+pH+sulphates+fixed.acidity+citric.acid+free.sulfur.dioxide2

################################Initial Model for Table S1######################
summary(model1)

################################Functional Residual for Initial Model###########

n<-nrow(model1@y)
probmodel1<-cbind.data.frame(rep(0,nrow(model1@y)),fitted(model1))
probrange1<-matrix(NA,nrow =nrow(model1@y) ,ncol=2)
y<-as.numeric(apply(model1@y, 1, function(t) colnames(model1@y)[which.max(t)]))
ordery<-y-min(y)+1
# range for the functional residual
for (i in 1:length(y)) {
  probrange1[i,]<-c(sum(probmodel1[i,1:ordery[i]]),sum(probmodel1[i,1:(ordery[i]+1)]))
}
#######heatmap part

prenumber<-matrix(NA,nrow=n,ncol = 11)
for (h in 1:n) {
  for (a in 1:ncol(prenumber)) {
    prenumber[h,a]<-probrange1[h,1]+(probrange1[h,2]-probrange1[h,1])/10*(a-1)
  }
}
prenumber_vector<-as.vector(prenumber)


q_prenumber_vector<-qnorm(prenumber_vector)

repchlorides<-rep(whitewine$chlorides,11)
repfixed.acidity<-rep(whitewine$fixed.acidity,11)
repvolatile.acidity<-rep(whitewine$volatile.acidity,11)
repcitric.acid<-rep(whitewine$citric.acid,11)
represidual.sugar<-rep(whitewine$residual.sugar,11)
repfree.sulfur.dioxide<-rep(whitewine$free.sulfur.dioxide,11)
reptotal.sulfur.dioxide<-rep(whitewine$total.sulfur.dioxide,11)
repdensity<-rep(whitewine$density,11)
reppH<-rep(whitewine$pH,11)
repsulphates<-rep(whitewine$sulphates,11)
repalcohol<-rep(whitewine$alcohol,11)



heatmapdata<-cbind.data.frame(repchlorides,repfixed.acidity,repvolatile.acidity,
                              repcitric.acid,represidual.sugar,
                              repfree.sulfur.dioxide,reptotal.sulfur.dioxide,
                              repdensity,reppH,repalcohol,repsulphates,q_prenumber_vector)

heatmap1_norm<-ggplot(heatmapdata, aes(repchlorides,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("chlorides")+ylab("")+
  labs(title = "Functional residuals")

heatmap2_norm<-ggplot(heatmapdata, aes(repfixed.acidity,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(a) fixed.acidity")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap3_norm<-ggplot(heatmapdata, aes(repvolatile.acidity,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = " volatile.acidity")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap4_norm<-ggplot(heatmapdata, aes(repcitric.acid,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("citric.acid")+ylab("")+
  labs(title = "Functional residuals")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap5_norm<-ggplot(heatmapdata, aes(represidual.sugar,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(b) residual.sugar")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap6_norm<-ggplot(heatmapdata, aes(repfree.sulfur.dioxide,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-5,3)+
  labs(title = "free.sulfur.dioxide")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))



heatmap7_norm<-ggplot(heatmapdata, aes(reptotal.sulfur.dioxide,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("total.sulfur.dioxide")+ylab("")+
  labs(title = "Functional residuals")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap8_norm<-ggplot(heatmapdata, aes(repdensity,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(c) density")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap9_norm<-ggplot(heatmapdata, aes(reppH,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "pH")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap10_norm<-ggplot(heatmapdata, aes(repalcohol,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "alcohol")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap11_norm<-ggplot(heatmapdata, aes(repsulphates,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "sulphates")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

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
heatmapdata<-heatmapdata %>%
  filter(represidual.sugar<=50)

heatmapdata<-heatmapdata %>%
  filter(repdensity<1.01)

heatmapdata<-heatmapdata %>%
  filter(repfixed.acidity<11)

heatmap2_norm_v2<-ggplot(heatmapdata, aes(repfixed.acidity,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(a*) fixed.acidity")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap5_norm_v2<-ggplot(heatmapdata, aes(represidual.sugar,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(b*) residual.sugar")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap8_norm_v2<-ggplot(heatmapdata, aes(repdensity,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(c*) density")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
####################################################
#####################Figure S10#####################
####################################################

grid.arrange(heatmap2_norm,heatmap5_norm,heatmap8_norm,
             heatmap2_norm_v2,heatmap5_norm_v2,heatmap8_norm_v2,nrow=2)



#################### Add square term #######################
whitewine<-read.csv("../Real-data/winequality-white.csv",sep = ";")
whitewine<-whitewine %>%
  filter(residual.sugar<=50)

whitewine<-whitewine %>%
  filter(density<1.01)

whitewine<-whitewine %>%
  filter(fixed.acidity<11)
whitewine$free.sulfur.dioxide2<-whitewine$free.sulfur.dioxide^2

model2<- vglm(quality~volatile.acidity+
                alcohol+sulphates+fixed.acidity+
                residual.sugar+free.sulfur.dioxide+
                pH+density+free.sulfur.dioxide2,
              family=acat(reverse=TRUE, parallel=TRUE),data =whitewine)

#alcohol+volatile.acidity+residual.sugar+free.sulfur.dioxide+density+pH+sulphates+fixed.acidity+citric.acid+free.sulfur.dioxide2
################################Final Model for Table S1######################

summary(model2)

##try one by one AIC=11023.12
#chlorides+volatile.acidity+total.sulfur.dioxide+alcohol
#sulphates+fixed.acidity+citric.acid+residual.sugar+free.sulfur.dioxide+pH+density



n<-nrow(model2@y)
probmodel2<-cbind.data.frame(rep(0,nrow(model2@y)),fitted(model2))
probrange2<-matrix(NA,nrow =nrow(model2@y) ,ncol=2)
y<-as.numeric(apply(model2@y, 1, function(t) colnames(model2@y)[which.max(t)]))
ordery<-y-min(y)+1
for (i in 1:length(y)) {
  probrange2[i,]<-c(sum(probmodel2[i,1:ordery[i]]),sum(probmodel2[i,1:(ordery[i]+1)]))
}

prenumber2<-matrix(NA,nrow=n,ncol = 11)
for (h in 1:n) {
  for (a in 1:ncol(prenumber2)) {
    prenumber2[h,a]<-probrange2[h,1]+(probrange2[h,2]-probrange2[h,1])/10*(a-1)
  }
}
prenumber_vector2<-as.vector(prenumber2)
q_prenumber_vector2<-qnorm(prenumber_vector2)

repchlorides<-rep(whitewine$chlorides,11)
repfixed.acidity<-rep(whitewine$fixed.acidity,11)
repvolatile.acidity<-rep(whitewine$volatile.acidity,11)
repcitric.acid<-rep(whitewine$citric.acid,11)
represidual.sugar<-rep(whitewine$residual.sugar,11)
repfree.sulfur.dioxide<-rep(whitewine$free.sulfur.dioxide,11)
reptotal.sulfur.dioxide<-rep(whitewine$total.sulfur.dioxide,11)
repdensity<-rep(whitewine$density,11)
reppH<-rep(whitewine$pH,11)
repsulphates<-rep(whitewine$sulphates,11)
repalcohol<-rep(whitewine$alcohol,11)



heatmapdata2<-cbind.data.frame(repchlorides,repfixed.acidity,repvolatile.acidity,
                              repcitric.acid,represidual.sugar,
                              repfree.sulfur.dioxide,reptotal.sulfur.dioxide,
                              repdensity,reppH,repalcohol,repsulphates,q_prenumber_vector2)

heatmapbefore<-ggplot(heatmapdata, aes(repfree.sulfur.dioxide,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-5,5)+
  labs(title = "(a) before")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

heatmap_norm_after_quard<-ggplot(heatmapdata2, aes(repfree.sulfur.dioxide,q_prenumber_vector2)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-3,3)+
  labs(title = "(b) after")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

####################################################
#####################Figure S12#####################
####################################################
multiplot(heatmapbefore,heatmap_norm_after_quard,cols = 2) 

