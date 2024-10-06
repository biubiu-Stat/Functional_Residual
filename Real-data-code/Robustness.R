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


set.seed(6)


whitewine<-read.csv("../Real-data/winequality-white.csv",sep = ";")
subsample<-sample(c(1:nrow(whitewine)),size = 250,replace = F)
whitewine<-whitewine[subsample,]


#alcohol+volatile.acidity+residual.sugar+free.sulfur.dioxide+density+pH+sulphates+fixed.acidity+citric.acid
model1<- vglm(quality~volatile.acidity+
                alcohol+sulphates+fixed.acidity+
                residual.sugar+free.sulfur.dioxide+
                pH+density,
              family=acat(reverse=TRUE, parallel=TRUE),data =whitewine)

n<-nrow(model1@y)
probmodel1<-cbind.data.frame(rep(0,nrow(model1@y)),fitted(model1))
probrange1<-matrix(NA,nrow =nrow(model1@y) ,ncol=2)
y<-as.numeric(apply(model1@y, 1, function(t) colnames(model1@y)[which.max(t)]))
ordery<-y-min(y)+1
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
heatmap6_norm<-ggplot(heatmapdata, aes(repfree.sulfur.dioxide,q_prenumber_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-3,3)+
  labs(title = "(a) free.sulfur.dioxide")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


bikedata<-read.csv("../Real-data/hour.csv")


bikedata<-bikedata %>% 
  filter(yr==1)

bikedata<-bikedata[,!grepl("casual",colnames(bikedata))]
bikedata<-bikedata[,!grepl("registered",colnames(bikedata))]
bikedata<-bikedata[,!grepl("atemp",colnames(bikedata))]
bikedata<-bikedata[,!grepl("yr",colnames(bikedata))]
bikedata<-bikedata[,!grepl("dteday",colnames(bikedata))]
bikedata<-bikedata[,!grepl("holiday",colnames(bikedata))]
bikedata<-bikedata[,!grepl("instant",colnames(bikedata))]
bikedata<-bikedata[,!grepl("mnth",colnames(bikedata))]
bikedata<-bikedata[,!grepl("weekday",colnames(bikedata))]

bikedata$winter<-ifelse(bikedata$season==1,1,0)
bikedata<-bikedata[,!grepl("season",colnames(bikedata))]
subsample<-sample(c(1:nrow(bikedata)),size = 250,replace = F)
bikedata<-bikedata[subsample,]
model1<-glm(cnt~.,family = "poisson",data = bikedata)
fittedy1<-model1$fitted.values
n<-length(fittedy1)

######heatmap
bikerange<-matrix(NA,nrow=n,ncol = 2)
for (i in 1:n) {
  bikerange[i,]<-c(ppois(bikedata$cnt[i]-1,fittedy1[i]),
                   ppois(bikedata$cnt[i],fittedy1[i]))
}
pre<-matrix(NA,nrow=n,ncol = 11)
for (h in 1:n) {
  for (a in 1:ncol(pre)) {
    pre[h,a]<-bikerange[h,1]+(bikerange[h,2]-bikerange[h,1])/10*(a-1)
  }
}

qnumbervector<-qnorm(as.vector(pre))

repwinter<-rep(bikedata$winter,11)+runif(length(qnumbervector), -0.1, 0.1)
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

bikeheatmapdata$repwinter<-bikeheatmapdata$repwinter+runif(nrow(bikeheatmapdata), -0.05, 0.05)

heatmap2_norm<-ggplot(bikeheatmapdata, aes(rephr,qnumbervector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-20,20)+
  labs(title = "(b) hour")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))



##########additive model


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

###########

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
  labs(title = "(c) Fn-Fn plot")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

#####quasipoisson


model_gam_quasi<-gam(cnt~winter+s(hr)+workingday+weathersit+
                       s(temp)+s(hum)+s(windspeed),
                     family = quasipoisson,
                     data = bikedata)
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


heatmap2_norm_gam_quasi<-ggplot(gam_quasi_heatmapdata, aes(rephr,gam_quasi_vector)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(d) hour")+ylim(-3,3)+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


################################################################################
###################################Figure S24###################################
################################################################################
grid.arrange(heatmap6_norm,heatmap2_norm,
             tpoints1_gam,heatmap2_norm_gam_quasi,nrow=2)
