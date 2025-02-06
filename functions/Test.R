
##################################################################################
################example 1 & 2 Missing of the quadratic term#######################
##################################################################################
set.seed(3)
library(VGAM)
library(brglm2)
library(ggplot2)
library(gridExtra)


n<-1000
x1<-rnorm(1000,0,1)
x2<-x1^2

linearp1<-1.5*x1-x2# link only for link function. linear predictor
linearp2<-1.5+1.5*x1-x2
linearp3<--1+1.5*x1-x2
linearp4<-1+1.5*x1-x2

p1<-1/(1+exp(-linearp1)+exp(-(linearp1+linearp2))+exp(-(linearp1+linearp2+linearp3))+exp(-(linearp1+linearp2+linearp3+linearp4)))
p2<-exp(-linearp1)/(1+exp(-linearp1)+exp(-(linearp1+linearp2))+exp(-(linearp1+linearp2+linearp3))+exp(-(linearp1+linearp2+linearp3+linearp4)))
p3<-exp(-(linearp1+linearp2))/(1+exp(-linearp1)+exp(-(linearp1+linearp2))+exp(-(linearp1+linearp2+linearp3))+exp(-(linearp1+linearp2+linearp3+linearp4)))
p4<-exp(-(linearp1+linearp2+linearp3))/(1+exp(-linearp1)+exp(-(linearp1+linearp2))+exp(-(linearp1+linearp2+linearp3))+exp(-(linearp1+linearp2+linearp3+linearp4)))
p5<-exp(-(linearp1+linearp2+linearp3+linearp4))/(1+exp(-linearp1)+exp(-(linearp1+linearp2))+exp(-(linearp1+linearp2+linearp3))+exp(-(linearp1+linearp2+linearp3+linearp4)))
pr<-cbind(p1,p2,p3,p4,p5)

y<-c()
for (i in 1:length(x1)) {
  y[i] <- sample(c(1,2,3,4,5), 1, replace=TRUE, prob=pr[i,]) 
}
testdata<-cbind.data.frame(x1,x2,y)

model1<- vglm(testdata$y~testdata$x1,
              family=acat(reverse=TRUE, parallel=TRUE))

model2<- vglm(testdata$y~testdata$x1 + testdata$x2,
              family=acat(reverse=TRUE, parallel=TRUE))#correct model

range1fun<-fresiduals(model1)
range2fun<-fresiduals(model2)


p1_unif <- fresplot(range1,x1,scale="uniform")
p2_unif <- fresplot(range2,x1,scale="uniform")
p1_norm <- fresplot(range1,x1,scale="normal",xl=-2,xp=2)
p2_norm <- fresplot(range2,x1,scale="normal")

grid.arrange(p1_unif,p1_norm,p2_unif,p2_norm,ncol=2)


numbers.t <- matrix(NA, nrow = 1000, ncol = 101)

# Fill the matrix with the computed values based on Functional_residual and x
for (h in 1:1000) {
  for (a in 1:ncol(numbers.t)) {
    numbers.t[h, a] <- range1fun[h, 1] + (range1fun[h, 2] - range1fun[h, 1]) / 100 * (a - 1)
  }
}

# Replicate x to match the number of rows in numbers
x_101t <- rep(x1, 101)
mean(qnumbers_vt==qnumbers1v)
# Convert the matrix to a vector and then apply the normal quantile transformation
numbers_vt <- as.vector(numbers.t)
qnumbers_vt <- qnorm(numbers_vt)

# Combine the x11 and transformed values into a data frame
qnumberst <- cbind.data.frame(x_101t, qnumbers_vt)
numbers.t <- cbind.data.frame(x_101t, numbers_vt)

numbers.t==numbers1v


p1_normt<-ggplot(qnumberst, aes(x101t,qnumbers_vt)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density")+
  geom_hline(yintercept=0,linetype="dashed", color = "red")+
  xlab("X")+ylab("")+xlim(-2,2)+
  geom_smooth(method = "loess",se=FALSE)+
  labs(title = "(b) Functional residuals on the normal scale")+
  theme(plot.title = element_text(size=12),axis.title=element_text(size=12))

