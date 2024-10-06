
###################Count Data Examples##############################

##################################################################################
################Examples 5 and S1 (a) #######################################################
##################################################################################
library(Rmisc)
library(plyr)
library(ggplot2)
library(VGAM)
library(pscl)
library(ggpointdensity)
library(MASS)

set.seed(3)
n<-1000
px1<-rnorm(1000,0,1)
px2<-px1^2
plinearp<-1+0.2*px1+0.15*px2# link only for link function. linear predictor
plambda<-exp(plinearp)
py<-c()
for (i in 1:length(px1)) {
  py[i]<-rpois(1,plambda[i])
  
}
pdata<-cbind.data.frame(px1,px2,py)
summary(py)
pmodel1<-glm(py~px1,family = "poisson",data = pdata)
pmodel2<-glm(py~px1+px2,family = "poisson",data = pdata)
summary(pmodel1)

summary(pmodel2)
fittedy1<-pmodel1$fitted.values
fittedy2<-pmodel2$fitted.values
pjitter1<-c()
range1<-matrix(NA,nrow=1000,ncol = 2)
range2<-matrix(NA,nrow=1000,ncol = 2)

for (i in 1:length(px1)) {
  range1[i,]<-c(ppois(py[i]-1,fittedy1[i]),ppois(py[i],fittedy1[i]))
}

for (i in 1:length(px1)) {
  range2[i,]<-c(ppois(py[i]-1,fittedy2[i]),ppois(py[i],fittedy2[i]))
}
numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}
px11<-rep(px1,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers1<-cbind.data.frame(px11,qnumbers1v)
qnumbers2<-cbind.data.frame(px11,qnumbers2v)
numbers1<-cbind.data.frame(px11,numbers1v)
numbers2<-cbind.data.frame(px11,numbers2v)
deviance1<-resid(pmodel1,type="deviance")
pearson1<-resid(pmodel1,type="pearson")
deviance2<-resid(pmodel2,type="deviance")
pearson2<-resid(pmodel2,type="pearson")
d1<-cbind.data.frame(px1,deviance1)
d2<-cbind.data.frame(px1,deviance2)
p1<-cbind.data.frame(px1,pearson1)
p2<-cbind.data.frame(px1,pearson2)

p1_unif<-ggplot(numbers1, aes(px11,numbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+ylim(0,1)+
  xlab("X")+ylab("")+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  labs(title = "Functional residuals on the uniform scale")

p2_unif<-ggplot(numbers2, aes(px11,numbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  ylim(0,1)+
  xlab("X")+ylab("")+labs(title ="(a) Functional residuals on the uniform scale")
p2_norm<-ggplot(qnumbers2, aes(px11,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+ylim(-3,3)+
  xlab("X")+ylab("")+
  labs(title = "(b) Functional residuals on the normal scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p1_norm_quard<-ggplot(qnumbers1, aes(px11,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("X")+ylab("")+ylim(-3,3)+
  geom_smooth(method = "loess",se=FALSE)+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  labs(title = "(a)")

p1_deviance<-ggplot(d1, aes(x=px1, y=deviance1)) + 
  geom_point()+
  geom_smooth(method=loess, se=FALSE)+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  xlab("X")+ylab("")+labs(title="Deviance residuals")+ylim(-4,4)
p2_deviance<-ggplot(d2, aes(x=px1, y=deviance2)) + 
  geom_point()+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  geom_smooth(method=loess, se=FALSE)+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("X")+ylab("")+labs(title="(c) Deviance Residuals")+ylim(-4,4)

p1_pearson<-ggplot(d1, aes(x=px1, y=pearson1)) + 
  geom_point()+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  geom_smooth(method=loess, se=FALSE)+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("X")+ylab("")+labs(title="Pearson residuals")+ylim(-4,4)
p2_pearson<-ggplot(d2, aes(x=px1, y=pearson2)) + 
  geom_point()+theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+
  geom_smooth(method=loess, se=FALSE)+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("X")+ylab("")+labs(title="(d) Pearson residuals")+ylim(-4,4)

########################################################
#######################Figure S5########################
########################################################
multiplot(p2_unif,p2_deviance,p2_norm,p2_pearson,cols=2)

########################################################
#######################Figure S6########################
########################################################

t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-as.data.frame(range1[1:100,])
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest[h,1],max=range1fortest[h,2])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(b)")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))+ylim(0,1)

t2<-seq(0,1,0.001)
res2<-c()
meanres2<-c()

range2fortest<-as.data.frame(range2)
for (i in 1:length(t2)) {
  for (h in 1:nrow(range2fortest)) {
    res2[h]<-punif(t2[i],min=range2fortest[h,1],max=range2fortest[h,2])
  }
  meanres2[i]<-mean(res2)
}

resdata2<-cbind.data.frame(t2,meanres2)

tpoints2<-ggplot(resdata2, aes(x=t2, y=meanres2)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(a)")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


t3<-seq(0,1,0.001)
res3<-c()
meanres3<-c()
lessthan0list<-which(px1<0)
range3fortest<-as.data.frame(range2)[lessthan0list,]
for (i in 1:length(t3)) {
  for (h in 1:nrow(range3fortest)) {
    res3[h]<-punif(t3[i],min=range3fortest[h,1],max=range3fortest[h,2])
  }
  meanres3[i]<-mean(res3)
}

resdata3<-cbind.data.frame(t3,meanres3)

tpoints3<-ggplot(resdata3, aes(x=t3, y=meanres3)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(b)")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

multiplot(tpoints2,tpoints3,cols=2)


##################################################################################
#####################Example S1 Missing of the higher order term#################
##################################################################################
set.seed(3)
px1<-rnorm(1000,0,0.5)
n <- 1000
px2<-px1^2
px3<-px1^3

plinearp<-0.8-0.2*px1+0.5*px2-0.5*px3# link only for link function. linear predictor

plambda<-exp(plinearp)

py<-c()
for (i in 1:length(px1)) {
  py[i]<-rpois(1,plambda[i])
}

summary(py)
pmodel1<-glm(py~px1+px2,family = "poisson")
pmodel2<-glm(py~px1+px2+px3,family = "poisson")
fittedy1<-pmodel1$fitted.values
fittedy2<-pmodel2$fitted.values
pjitter1<-c()
range1<-matrix(NA,nrow=1000,ncol = 2)
range2<-matrix(NA,nrow=1000,ncol = 2)

for (i in 1:length(px1)) {
  range1[i,]<-c(ppois(py[i]-1,fittedy1[i]),ppois(py[i],fittedy1[i]))
}

for (i in 1:length(px1)) {
  range2[i,]<-c(ppois(py[i]-1,fittedy2[i]),ppois(py[i],fittedy2[i]))
}
numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}
px11<-rep(px1,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers1<-cbind.data.frame(px11,qnumbers1v)
qnumbers2<-cbind.data.frame(px11,qnumbers2v)
numbers1<-cbind.data.frame(px11,numbers1v)
numbers2<-cbind.data.frame(px11,numbers2v)

p1_unif<-ggplot(numbers1, aes(px11,numbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+ylim(0,1)+
  xlab("X")+ylab("")+
  labs(title = "Functional residuals on the uniform scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p2_unif<-ggplot(numbers2, aes(px11,numbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  ylim(0,1)+
  xlab("X")+ylab("")+labs(title ="Functional residuals on the uniform scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
p2_norm<-ggplot(qnumbers2, aes(px11,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+xlim(-3,3)+
  xlab("X")+ylab("")+
  labs(title = "Functional residuals on the normal scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p1_norm_cubic<-ggplot(qnumbers1, aes(px11,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("X")+ylab("")+ylim(-3,3)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(b)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

############################################################################
###########Figure S17 (a)quadratic term is not included from Example 5######
###########Figure S17 (b)cubic term is not included ########################
############################################################################
multiplot(p1_norm_quard,p1_norm_cubic,cols=2)


########################        end        ######################################
#################################################################################




#################################################################################
############    Example 6    zero inflation        ##############################
#################################################################################
library(VGAM)
library(pscl)
library(ggplot2)
set.seed(3)
n<-1000
px1<-rnorm(1000,0,0.8)

plinearp<-1+1*px1# only for link function. linear predictor
plambda<-exp(plinearp)
p0<-exp(1+0.2*px1)/(exp(1+0.2*px1)+1)
py<-rzipois(n,lambda=plambda,pstr0 = p0)
pdata<-cbind.data.frame(px1,py)

pmodel1<-glm(py~px1,family = poisson(link=log),data = pdata)#Regular model

pmodel2<-zeroinfl(py~px1,data = pdata)# Zero-inflated Poisson model


fittedy1<-predict(pmodel1,type="response")

#Fitted value for Zero-inflated model
u2<-exp(coef(pmodel2)[1]+coef(pmodel2)[2]*px1)
pi2<-exp(coef(pmodel2)[3]+coef(pmodel2)[4]*px1)/(exp(coef(pmodel2)[3]+coef(pmodel2)[4]*px1)+1)

range1 <- as.data.frame(cbind(ppois(py - 1, fittedy1), ppois(py, fittedy1)))


range2 <- as.data.frame(cbind(pzipois(py - 1, lambda = u2,pstr0 = pi2)
                              , pzipois(py, lambda = u2,pstr0 = pi2)))



numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}
px11<-rep(px1,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers1<-cbind.data.frame(px11,qnumbers1v)
qnumbers2<-cbind.data.frame(px11,qnumbers2v)
numbers1<-cbind.data.frame(px11,numbers1v)
numbers2<-cbind.data.frame(px11,numbers2v)


p2_norm<-ggplot(qnumbers2, aes(px11,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+
  labs(title = "(b) Zero-Inflated Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p1_norm<-ggplot(qnumbers1, aes(px11,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("")+ylab("")+ylim(-2,2)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(a) Regular Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

## Fn-Fn Plots

t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-range1[1:1000,]
range2fortest<-range2[1:1000,]
range2fortest<-as.data.frame(range2fortest)
range1fortest<-as.data.frame(range1fortest)
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest$V1[h],max=range1fortest$V2[h])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(c) Regular Poisson model")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


res2<-c()
meanres2<-c()
for (i in 1:length(t1)) {
  for (h in 1:nrow(range2fortest)) {
    res2[h]<-punif(t1[i],min=range2fortest$V1[h],max=range2fortest$V2[h])
  }
  meanres2[i]<-mean(res2)
}

resdata2<-cbind.data.frame(t1,meanres2)

tpoints2<-ggplot(resdata2, aes(x=t1, y=meanres2)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(d) Zero-Inflated Poisson model" )+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

########################################################
#######################Figure S7########################
########################################################
multiplot(p1_norm,tpoints1,p2_norm,tpoints2,cols=2)







#################################################################
############  Example 7 Dispersed Poisson Model #################
#################################################################
library(ggplot2)
library(MASS)
library(gridExtra)
set.seed(3)
n<-1000
px1<-rnorm(1000,0,1)
rpois.od<-function (n, lambda,d) {
  if (d==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}
plinearp<-1.2+1.3*px1# link only for link function. linear predictor
f<-7
plambda<-exp(plinearp)
mean(plambda)
py<-rpois.od(n,plambda,f)
pdata<-cbind.data.frame(px1,py)
pmodel1<-glm(py~px1,family = poisson(link=log),data = pdata)
pmodel2<-glm(py~px1,family = quasipoisson)
s1<-summary(pmodel1)
s2<-summary(pmodel2)

s1$dispersion # default for regular Poisson model is 1
s2$dispersion

fittedy1<-predict(pmodel1,type="response")
fittedy2<-predict(pmodel2,type="response")

range1<-matrix(NA,nrow=1000,ncol = 2)
range2<-matrix(NA,nrow=1000,ncol = 2)
for (i in 1:length(px1)) {
  range1[i,]<-c(ppois(py[i]-1,fittedy1[i]),ppois(py[i],fittedy1[i]))
}

for (i in 1:length(px1)) {
  range2[i,]<-c(pnbinom(py[i]-1,size=fittedy2[i]/s2$dispersion,mu=fittedy2[i]),
                pnbinom(py[i],size=fittedy2[i]/s2$dispersion,mu=fittedy2[i]))
}

numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}
px11<-rep(px1,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers1<-cbind.data.frame(px11,qnumbers1v)
qnumbers2<-cbind.data.frame(px11,qnumbers2v)
numbers1<-cbind.data.frame(px11,numbers1v)
numbers2<-cbind.data.frame(px11,numbers2v)

p1_unif<-ggplot(numbers1, aes(px11,numbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+ylim(0,1)+
  xlab("X")+ylab("")+
  labs(title = "Functional residuals on the uniform scale")

p2_unif<-ggplot(numbers2, aes(px11,numbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  ylim(0,1)+
  xlab("X")+ylab("")+labs(title ="Functional residuals on the uniform scale")
p2_norm<-ggplot(qnumbers2, aes(px11,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-3,3)+
  labs(title = "(b) Quasi-Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p1_norm<-ggplot(qnumbers1, aes(px11,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("")+ylab("")+ylim(-3,3)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(a) Regular Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

# Fn-Fn Plot

t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-range1[1:1000,]
range2fortest<-range2[1:1000,]
range2fortest<-as.data.frame(range2fortest)
range1fortest<-as.data.frame(range1fortest)
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest$V1[h],max=range1fortest$V2[h])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(c) Regular Poisson model")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))



res2<-c()
meanres2<-c()
for (i in 1:length(t1)) {
  for (h in 1:nrow(range2fortest)) {
    res2[h]<-punif(t1[i],min=range2fortest$V1[h],max=range2fortest$V2[h])
  }
  meanres2[i]<-mean(res2)
}

resdata2<-cbind.data.frame(t1,meanres2)

tpoints2<-ggplot(resdata2, aes(x=t1, y=meanres2)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(d) Quasi-Poisson model" )+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
########################################################
#######################Figure S8########################
########################################################
grid.arrange(p1_norm,p2_norm,tpoints1,tpoints2,nrow=2)






################################################################# 
############# Example 8 Semi-parametric Poisson model############
################################################################# 
library(ggplot2)
library(gridExtra)
set.seed(3)
n<-1000
px1<-rnorm(1000,0,1)

rpois.od<-function (n, lambda,d) {
  if (d==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}
plinearp<-1.2+1.3*sin(px1)-0.8*px1# link only for link function. linear predictor
f<-7 # over dispersion parameter

plambda<-exp(plinearp)
mean(plambda)
py<-rpois.od(n,plambda,f)
pdata<-cbind.data.frame(px1,py)

pmodel1<-glm(py~px1,family = poisson(link=log),data = pdata)
# regular Poisson model
pmodel2<-gam(py~s(px1),family = poisson(link=log),data = pdata)
# Generalized additive Poisson model
pmodel3_gam<-gam(py~s(px1),family = quasipoisson,data = pdata)
# Generalized additive Quasi-Poisson model

s1<-summary(pmodel1)
s2<-summary(pmodel2)
s3<-summary(pmodel3_gam)


fittedy1<-predict(pmodel1,type="response")
fittedy2<-predict(pmodel2,type="response")
fittedy3<-predict(pmodel3_gam,type="response")

range1<-matrix(NA,nrow=1000,ncol = 2)
range2<-matrix(NA,nrow=1000,ncol = 2)
range3<-matrix(NA,nrow=1000,ncol = 2)

for (i in 1:length(px1)) {
  range1[i,]<-c(ppois(py[i]-1,fittedy1[i]),ppois(py[i],fittedy1[i]))
}

for (i in 1:length(px1)) {
  range2[i,]<-c(pnbinom(py[i]-1,size=fittedy2[i]/s2$dispersion,mu=fittedy2[i]),
                pnbinom(py[i],size=fittedy2[i]/s2$dispersion,mu=fittedy2[i]))
}

for (i in 1:length(px1)) {
  range3[i,]<-c(pnbinom(py[i]-1,size=fittedy3[i]/s3$dispersion,mu=fittedy3[i]),
                pnbinom(py[i],size=fittedy3[i]/s3$dispersion,mu=fittedy3[i]))
}

numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}

numbers3<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers3)) {
    numbers3[h,a]<-range3[h,1]+(range3[h,2]-range3[h,1])/100*(a-1)
  }
}

px11<-rep(px1,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
numbers3v<-as.vector(numbers3)

qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers3v<-qnorm(numbers3v)

qnumbers1<-cbind.data.frame(px11,qnumbers1v)
qnumbers2<-cbind.data.frame(px11,qnumbers2v)
qnumbers3<-cbind.data.frame(px11,qnumbers3v)

numbers1<-cbind.data.frame(px11,numbers1v)
numbers2<-cbind.data.frame(px11,numbers2v)

# Functional residual density plots

p2_norm<-ggplot(qnumbers2, aes(px11,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab("")+ylab("")+ylim(-3,3)+
  labs(title = "(b) Generalized additive Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p1_norm<-ggplot(qnumbers1, aes(px11,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("")+ylab("")+ylim(-3,3)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(a) Regular Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p3_norm<-ggplot(qnumbers3, aes(px11,qnumbers3v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("")+ylab("")+ylim(-3,3)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(c) Generalized additive quasi-Poisson model")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

### Fn-Fn Plots 

t1<-seq(0,1,0.001)
res1<-c()
meanres1<-c()

range1fortest<-range1[1:1000,]
range2fortest<-range2[1:1000,]
range3fortest<-range3[1:1000,]

range2fortest<-as.data.frame(range2fortest)
range1fortest<-as.data.frame(range1fortest)
range3fortest<-as.data.frame(range3fortest)
for (i in 1:length(t1)) {
  for (h in 1:nrow(range1fortest)) {
    res1[h]<-punif(t1[i],min=range1fortest$V1[h],max=range1fortest$V2[h])
  }
  meanres1[i]<-mean(res1)
}

resdata1<-cbind.data.frame(t1,meanres1)

tpoints1<-ggplot(resdata1, aes(x=t1, y=meanres1)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(d) Regular Poisson model")+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))



res2<-c()
meanres2<-c()
for (i in 1:length(t1)) {
  for (h in 1:nrow(range2fortest)) {
    res2[h]<-punif(t1[i],min=range2fortest$V1[h],max=range2fortest$V2[h])
  }
  meanres2[i]<-mean(res2)
}

resdata2<-cbind.data.frame(t1,meanres2)

tpoints2<-ggplot(resdata2, aes(x=t1, y=meanres2)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(e) Generalized additive Poisson model" )+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

res3<-c()
meanres3<-c()
for (i in 1:length(t1)) {
  for (h in 1:nrow(range3fortest)) {
    res3[h]<-punif(t1[i],min=range3fortest$V1[h],max=range3fortest$V2[h])
  }
  meanres3[i]<-mean(res3)
}

resdata3<-cbind.data.frame(t1,meanres3)

tpoints3<-ggplot(resdata3, aes(x=t1, y=meanres3)) + 
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype="dashed", color = "red")+
  labs(title = "(f) Generalized additive quasi-Poisson model" )+ylab(expression(bar(Res)(t)))+xlab("t")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

########################################################
#######################Figure S9########################
########################################################
grid.arrange(p1_norm,p2_norm,p3_norm,tpoints1,tpoints2,tpoints3,nrow=2)









##################################################################################
######################Example S2 Missing of a covariate############################
##################################################################################
set.seed(3)
n<-1000
px1<-rnorm(1000,0,0.8)

px2<-rnorm(1000,-1,1)
px3<-rnorm(1000,0.8,0.9)
plinearp<-0.5+0.25*px1+0.5*px2# link only for link function. linear predictor
plambda<-exp(plinearp)

py<-c()
for (i in 1:length(px1)) {
  py[i]<-rpois(1,plambda[i])
}


pmodel1<-glm(py~px1,family = "poisson")
pmodel2<-glm(py~px1+px2,family = "poisson")
fittedy1<-pmodel1$fitted.values
fittedy2<-pmodel2$fitted.values
pjitter1<-c()
range1<-matrix(NA,nrow=1000,ncol = 2)
range2<-matrix(NA,nrow=1000,ncol = 2)

for (i in 1:length(px1)) {
  range1[i,]<-c(ppois(py[i]-1,fittedy1[i]),ppois(py[i],fittedy1[i]))
}

for (i in 1:length(px1)) {
  range2[i,]<-c(ppois(py[i]-1,fittedy2[i]),ppois(py[i],fittedy2[i]))
}
numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}
px11<-rep(px1,101)
px22<-rep(px2,101)
px33<-rep(px3,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers1<-cbind.data.frame(px22,qnumbers1v)
qnumbers2<-cbind.data.frame(px33,qnumbers1v)
numbers1<-cbind.data.frame(px22,numbers1v)
numbers2<-cbind.data.frame(px33,numbers1v)

p1_unif<-ggplot(numbers1, aes(px22,numbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+ylim(0,1)+
  xlab(expression(X))+ylab("")+
  labs(title = "Functional residuals on the uniform scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p2_unif<-ggplot(numbers2, aes(px33,numbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  ylim(0,1)+
  xlab(expression(X))+ylab("")+labs(title ="Functional residuals on the uniform scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
p2_norm<-ggplot(qnumbers2, aes(px33,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  xlab(expression(X[3]))+ylab("")+
  labs(title = "(b)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p1_norm<-ggplot(qnumbers1, aes(px22,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab(expression(X[2]))+ylab("")+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(a)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

##########################################################
##################### Figure S18 #########################
##########################################################
multiplot(p1_norm,p2_norm,cols = 2)


#################################################################################
########################        end        ######################################
#################################################################################



##################################################################################
######################Example S3 Missing of the interaction term##################
##################################################################################
set.seed(33)
n<-1000
px1<-rnorm(n,0.5,1)

px2<-rnorm(n,-1,0.7)
px3<-px1*px2
plinearp<--0.1+0.8*px1-0.5*px2+0.6*px3# link only for link function. linear predictor
plambda<-exp(plinearp)
summary(plambda)
py<-c()
for (i in 1:length(px1)) {
  py[i]<-rpois(1,plambda[i])
}

pmodel1<-glm(py~px1+px2,family = "poisson")
pmodel2<-glm(py~px1+px2+px3,family = "poisson")
fittedy1<-pmodel1$fitted.values
fittedy2<-pmodel2$fitted.values
pjitter1<-c()
range1<-matrix(NA,nrow=1000,ncol = 2)
range2<-matrix(NA,nrow=1000,ncol = 2)

for (i in 1:length(px1)) {
  range1[i,]<-c(ppois(py[i]-1,fittedy1[i]),ppois(py[i],fittedy1[i]))
}

for (i in 1:length(px1)) {
  range2[i,]<-c(ppois(py[i]-1,fittedy2[i]),ppois(py[i],fittedy2[i]))
}
numbers1<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers1)) {
    numbers1[h,a]<-range1[h,1]+(range1[h,2]-range1[h,1])/100*(a-1)
  }
}
numbers2<-matrix(NA,nrow=n,ncol = 101)
for (h in 1:n) {
  for (a in 1:ncol(numbers2)) {
    numbers2[h,a]<-range2[h,1]+(range2[h,2]-range2[h,1])/100*(a-1)
  }
}

px33<-rep(px3,101)
numbers1v<-as.vector(numbers1)
numbers2v<-as.vector(numbers2)
qnumbers1v<-qnorm(numbers1v)
qnumbers2v<-qnorm(numbers2v)
qnumbers1<-cbind.data.frame(px33,qnumbers1v)
qnumbers2<-cbind.data.frame(px33,qnumbers2v)
numbers1<-cbind.data.frame(px33,numbers1v)
numbers2<-cbind.data.frame(px33,numbers2v)

p1_unif<-ggplot(numbers1, aes(px33,numbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+ylim(0,1)+
  xlab(expression(X[1]*X[2]))+ylab("")+
  labs(title = "Functional residuals on the uniform scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

p2_unif<-ggplot(numbers2, aes(px33,numbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0.5,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+
  ylim(0,1)+
  xlab(expression(X[1]*X[2]))+ylab("")+labs(title ="Functional residuals on the uniform scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
p2_norm<-ggplot(qnumbers2, aes(px33,qnumbers2v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  geom_smooth(method = "loess",se=FALSE)+xlim(-4,4)+
  xlab(expression(X[1]*X[2]))+ylab("")+
  labs(title = "(b)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))
qnumbers1<-cbind.data.frame(px33,qnumbers1v)
p1_norm<-ggplot(qnumbers1, aes(px33,qnumbers1v)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab(expression(X[1]*X[2]))+ylab("")+xlim(-4,4)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(a)")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))


##########################################################
##################### Figure S19 #########################
##########################################################
multiplot(p1_norm,p2_norm,cols=2)



#################################################################################
########################        end        ######################################
#################################################################################




