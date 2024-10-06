##########FnFnPlot with interval
library(MASS)
library(ggplot2)
library(gridExtra)
library(VGAM)
library(pscl)
library(ggplot2)

set.seed(3)

########## Simulation for Overdispersion ##############

rpois.od<-function (n, lambda,d) {
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}


n<-nn<-1000
px1<-rnorm(1000,0,1)
t1<-seq(0,1,0.02)
meanres1B <- matrix(NA, nrow = 2000, ncol = 51)
meanres2B <- matrix(NA, nrow = 2000, ncol = 51)



plinearp<-1.2+1.3*px1# for link function. linear predictor
f<-7
plambda<-exp(plinearp)
mean(plambda)
py<-rpois.od(n,plambda,f)
pdata<-cbind.data.frame(px1,py)
pmodel1<-glm(py~px1,family = poisson(link=log),data = pdata)
pmodel2<-glm(py~px1,family = quasipoisson)
s1<-summary(pmodel1)
s2<-summary(pmodel2)

s1$dispersion
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

range1<-as.data.frame(range1)
range2<-as.data.frame(range2)

meanres1 <- sapply(t1, function(t) {
  res1 <- punif(t, min = range1$V1, max = range1$V2)
  mean(res1)
})

meanres2 <- sapply(t1, function(t) {
  res2 <- punif(t, min = range2$V1, max = range2$V2)
  mean(res2)
})

########## Bootstrapping for probable envelop ##############

for (B in 1:2000) {
  index<-sample(1:n,nn,replace = T)
  pdataB<-pdata[index,]
  model1B<-glm(py~px1,family = "poisson",data = pdataB)
  model2B<-glm(py~px1,family = quasipoisson,data = pdataB)
  fittedy1B<-model1B$fitted.values
  fittedy2B<-model2B$fitted.values
  s2B<-summary(model2B)
  range1B <- as.data.frame(cbind(ppois(pdataB$py - 1, fittedy1B), 
                                 ppois(pdataB$py, fittedy1B)))
  range2B<-matrix(NA,nrow=n,ncol = 2)
  for (i in 1:nn) {
    range2B[i,]<-c(pnbinom(pdataB$py[i]-1,size=fittedy2B[i]/s2B$dispersion,
                           mu=fittedy2B[i]),
                   pnbinom(pdataB$py[i],size=fittedy2B[i]/s2B$dispersion,
                           mu=fittedy2B[i]))
  }
  range2B<-as.data.frame(range2B)
  
  
  meanres1B[B,] <- sapply(t1, function(t) {
    res1B <- punif(t, min = range1B$V1, max = range1B$V2)
    mean(res1B)
  })
  
  meanres1B[B,]<-meanres1B[B,]-t1
  
  meanres2B[B,] <- sapply(t1, function(t) {
    res2B <- punif(t, min = range2B$V1, max = range2B$V2)
    mean(res2B)
  })
  meanres2B[B,]<-meanres2B[B,]-t1
}
meanres1B<-as.data.frame(meanres1B)
meanres2B<-as.data.frame(meanres2B)


interval1<-matrix(NA,nrow=51,ncol=2)
interval2<-matrix(NA,nrow=51,ncol=2)
for (j in 1:51) {
  interval1[j,]<-quantile(meanres1B[,j],probs = c(0.025,0.975))
  interval2[j,]<-quantile(meanres2B[,j],probs = c(0.025,0.975))
}

interval1<-as.data.frame(interval1)
interval2<-as.data.frame(interval2)

interval1<-cbind.data.frame(t1,interval1,meanres1)
interval2<-cbind.data.frame(t1,interval2,meanres2)



intervalLine1 <- ggplot(interval1, aes(t1)) +
  geom_line(aes(y = V1 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = V2 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = meanres1), colour = "black", size = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "(a) Regular Poisson model",
       x = "t")+ylab(expression(bar(Res)(t))) 


intervalLine2<-ggplot(interval2, aes(t1)) +
  geom_line(aes(y = V1 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = V2 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = meanres2), colour = "black", size = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "(b) Quasi-Poisson model",
       x = "t")+ylab(expression(bar(Res)(t)))




########### Simulation for Zero-inflated #############

set.seed(3)
n<-nn<-1000
t1<-seq(0,1,0.02)
px1<-rnorm(1000,0,0.8)

meanres1Bz <- matrix(NA, nrow = 2000, ncol = 51)
meanres2Bz <- matrix(NA, nrow = 2000, ncol = 51)



plinearp<-1+1*px1# link only for link function. linear predictor
plambda<-exp(plinearp)
p0<-exp(1+0.2*px1)/(exp(1+0.2*px1)+1)
py<-rzipois(n,lambda=plambda,pstr0 = p0)
pdata<-cbind.data.frame(px1,py)

pmodel1<-glm(py~px1,family = poisson(link=log),data = pdata)
pmodel2<-zeroinfl(py~px1,data = pdata)

fittedy1<-predict(pmodel1,type="response")
u2<-exp(coef(pmodel2)[1]+coef(pmodel2)[2]*px1)
pi2<-exp(coef(pmodel2)[3]+coef(pmodel2)[4]*px1)/(exp(coef(pmodel2)[3]+coef(pmodel2)[4]*px1)+1)


range1z <- as.data.frame(cbind(ppois(py - 1, fittedy1), ppois(py, fittedy1)))


range2z <- as.data.frame(cbind(pzipois(py - 1, lambda = u2,pstr0 = pi2)
                              , pzipois(py, lambda = u2,pstr0 = pi2)))

meanres1z <- sapply(t1, function(t) {
  res1z <- punif(t, min = range1z$V1, max = range1z$V2)
  mean(res1z)
})

meanres2z <- sapply(t1, function(t) {
  res2z <- punif(t, min = range2z$V1, max = range2z$V2)
  mean(res2z)
})
########## Bootstrapping for probable envelop ##############

for (B in 1:2000) {
  indexz<-sample(1:n,nn,replace = T)
  pdataBz<-pdata[indexz,]
  model1B<-glm(py~px1,family = "poisson",data = pdataBz)
  model2B<-zeroinfl(py~px1,data = pdataBz)
  
  fittedy1B<-model1B$fitted.values
  
  
  range1Bz <- as.data.frame(cbind(ppois(pdataBz$py - 1, fittedy1B), 
                                 ppois(pdataBz$py, fittedy1B)))
  u2B<-exp(coef(model2B)[1]+coef(model2B)[2]*pdataBz$px1)
  pi2B<-exp(coef(model2B)[3]+
              coef(model2B)[4]*pdataBz$px1)/(exp(coef(model2B)[3]+
                                                  coef(model2B)[4]*pdataBz$px1)+1)
  
  range2Bz <- as.data.frame(cbind(pzipois(pdataBz$py - 1, lambda = u2B,pstr0 = pi2B)
                                 , pzipois(pdataBz$py, lambda = u2B,pstr0 = pi2B)))
  
  
  
  meanres1Bz[B,] <- sapply(t1, function(t) {
    res1Bz <- punif(t, min = range1Bz$V1, max = range1Bz$V2)
    mean(res1Bz)
  })
  
  meanres1Bz[B,]<-meanres1Bz[B,]-t1
  
  meanres2Bz[B,] <- sapply(t1, function(t) {
    res2Bz <- punif(t, min = range2Bz$V1, max = range2Bz$V2)
    mean(res2Bz)
  })
  meanres2Bz[B,]<-meanres2Bz[B,]-t1
}

meanres1Bz<-as.data.frame(meanres1Bz)
meanres2Bz<-as.data.frame(meanres2Bz)



interval1z<-matrix(NA,nrow=51,ncol=2)
interval2z<-matrix(NA,nrow=51,ncol=2)
for (j in 1:51) {
  interval1z[j,]<-quantile(meanres1Bz[,j],probs = c(0.025,0.975))
  interval2z[j,]<-quantile(meanres2Bz[,j],probs = c(0.025,0.975))
}

interval1z<-as.data.frame(interval1z)
interval2z<-as.data.frame(interval2z)

interval1z<-cbind.data.frame(t1,interval1z,meanres1z)
interval2z<-cbind.data.frame(t1,interval2z,meanres1z)



intervalLine3 <- ggplot(interval1z, aes(t1)) +
  geom_line(aes(y = V1 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = V2 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = meanres1z), colour = "black", size = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "(c) Regular Poisson model",
       x = "t") +ylab(expression(bar(Res)(t)))


intervalLine4<-ggplot(interval2z, aes(t1)) +
  geom_line(aes(y = V1 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = V2 + t1), colour = "black", linetype = "dotted") +
  geom_line(aes(y = meanres2z), colour = "black", size = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "(d) Zero-inflated Poisson model",
       x = "t") +ylab(expression(bar(Res)(t)))


###################### Figure S23 #############################################

grid.arrange(intervalLine1,intervalLine2,intervalLine3,intervalLine4,nrow=2)

