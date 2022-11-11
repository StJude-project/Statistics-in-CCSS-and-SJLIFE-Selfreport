

rm(list=ls())

set.seed(5)

library(plyr)
library(dplyr)
library(mvtnorm)
library(Rmpfr)
library(stringr)
library(geepack)
library(sandwich)


# Data Generation

NVcov<-read.csv("NVdata.csv")
Vcov<-read.csv("Vdata.csv")

pos <- match(Vcov$ID, NVcov$ID)

ndata<- NVcov[pos,]

#We generate samples of having CHC outcome from a Binomial distribution with success probabilities P(Y |X) = X?? and P(??,??)(S|X)  for clinically assessed and self-reported outcomes respectively, where P(??,??)(S|X) = ??? y=0,1 P?? (S|Y = y, X)P?? (Y = y|X) P (S|Y, X) = ??0 + ??1Y + Z??.
# Calculate the probabilities of y=1  (P(Y=1 |X) = exp(X.beta)/1+exp(X.beta))

# True parameters
Beta<- c(-7,-0.4,1.5,2,1.3,0.1)
th0<--6
th1<-3
Theta<-c(th0,th1,-0.4,-1,-1,1,1,2,1,0.01)

# covariate names for y model
cov.1 <- c( "hrtdose1","hrtdose2","hrtdose3", "hypten","age")

# covariate names for s model(excluding interactions)
cov.2 <- c("y" ,  "crtdose1", "crtdose2","crtdose3", "hrtdose1","hrtdose2","hrtdose3", "hypten","age")


#============ Generate Y =============
vdata1<-NVcov[pos,]
vdata2<-cbind(intercept=1, vdata1[, cov.1])
n.v<-dim(vdata2)[1]
vdata3<-as.matrix(vdata2)

 Beta<-as.numeric(Beta)
 Py<-1 /(1+exp(-(vdata3 %*% Beta)))
 
 
#========== Generate S=============      
# calculate P(S|Y,X) = exp(X.Theta)/1+exp(X.Theta))

n.T<-dim(NVcov)[1]
covdata<- NVcov[,c(  "crtdose1",   "crtdose2" ,  "crtdose3" , "hrtdose1" ,  "hrtdose2", "hrtdose3",  "hypten"  ,"age")]

T.vdata<-cbind(intercept=1, NVcov[, cov.1])
T.vdata1<-as.matrix(T.vdata)
Beta<-as.numeric(Beta)
Py<-1 /(1+exp(-(T.vdata1 %*% Beta)))
y1<-NA
# -generate y - Bernoulli distribution 
for (i in 1:n.T ){
  y1[i]<- rbinom(1, 1,Py[i])
}

tot.dat<-cbind(intercept=1, y1, covdata)

T.data<-as.matrix(tot.dat)
Ps<-exp(T.data %*% Theta) /(1+exp(T.data %*% Theta))

s<-NA
# -generate y - Bernoulli distribution 
for (i in 1:n.T ){
  s[i]<- rbinom(1, 1,Ps[i])
}

tot.dat$s<- s

ID<-NVcov$ID
tot.data<-cbind(ID,tot.dat)


# ====extract the s outcome from total sample for the validated sample
vdata <- tot.data[pos,]
colnames(vdata)[colnames(vdata)=="y1"] <- "y"
y<-vdata$y
nvdata1 <- tot.data[-pos,]
drop <- c("y1")
nvdata = nvdata1[,!(names(nvdata1) %in% drop)]

#======== Function to list all needed data ========================================
dat.prep <- function(vdat=vdata, nvdat=nvdata, cov1=cov.1, cov2=cov.2){
  
  
  # matrix for validated sample, model for y
  mv.y <- cbind(intercept=1, vdat[, cov1])
  
  
  mv.s <- cbind(intercept=1, vdat[, cov2])
  
  # matrix for non-validates sample, model for y
  mn.y <- cbind(intercept=1, nvdat[, cov1])
  
  # matrix for non-validates sample, model for s
  ### y=1
  nvdat.y1 <- nvdat
  nvdat.y1$y <- 1
 
  mn.s1 <- cbind(intercept=1, nvdat.y1[, cov2])
  
  
  ### y=0
  nvdat.y0 <- nvdat
  nvdat.y0$y <- 0
  
  mn.s2 <- cbind(intercept=1, nvdat.y0[, cov2])
  
  
  # combine 5 matrices
  odata <- list( cov1=cov1, cov2=cov2,
                 y=y, s=vdat$s, s2=nvdat$s, 
                 mv.y=mv.y,   # matrix for validated, model for y
                 mv.s=mv.s,   # matrix for validated, model for s
                 mn.y=mn.y,   # matrix for non-validated, model for y
                 mn.s1=mn.s1, # matrix for non-validated, model for s, when y=1
                 mn.s2=mn.s2  # matrix for non-validated, model for s, when y=0
  )
  
}

odata <- dat.prep()



### $\color{red}{\text{Misclassification Table for the overlap sub-sample}}$


# data for double-sampling
dat.ds <- left_join(vdata[,c("ID", "s")], vdata[,c("ID", "y")])
tb <- xtabs(~ y+s, dat.ds)
tb <- matrix(c(tb[2,2], tb[1,2], tb[2,1], tb[1,1]), nrow=2) # change the order so it can be used in double.sample()

tb <- as.data.frame(tb)
colnames(tb) <- c("s=1", "s=0")
rownames(tb) <- c("y=1", "y=0")

  tb

  V_Y.NV_N<- tb[1,2]
V_N.NV_Y<- tb[2,1]

print(paste("# of Yes in Validated sample, No in non-validated sample   ",  V_Y.NV_N))
print(paste("# of No in Validated sample, Yes in non-validated sample   ", V_N.NV_Y))


####============Pepe's method=========================########
  
# log likelihood function

loglik.fnPepe <- function(par, data=odata){
  # parameters
  B <- par[1: ncol(data$mv.y)]
  Th <- par[(ncol(data$mv.y)+1): length(par)]
  y <- unlist(data$y)
  s <- unlist(data$s)
  s2 <- unlist(data$s2)
  mv.y <- matrix(unlist(data$mv.y), ncol=ncol(data$mv.y), byrow=F)
  mv.s <- matrix(unlist(data$mv.s), ncol=ncol(data$mv.s), byrow=F)
  mn.y <- matrix(unlist(data$mn.y), ncol=ncol(data$mn.y), byrow=F)
  mn.s1 <- matrix(unlist(data$mn.s1), ncol=ncol(data$mn.s1), byrow=F)
  mn.s2 <- matrix(unlist(data$mn.s2), ncol=ncol(data$mn.s1), byrow=F)
  
  ## lglik1 
  l1 <- rep(NA, length(y))
  
  eta1 <- mv.y %*% B
  
  lglik1 <- sum(eta1*y-log(1+exp(eta1)))
  
  
  ## lglik2
  l2  <- rep(NA, length(y))
  
  eta2 <- mv.s %*% Th
  
  lglik2 <-sum(eta2*s-log(1+exp(eta2)))
  
  ## lglik3
  l3 <-   rep(NA, length(s2)) 
  
  eta3 <- mn.y %*% B # eta3 has same form as eta1
  
  eta2_y1 <- mn.s1 %*% Th #y=1
  
  eta2_y0 <- mn.s2 %*% Th #y=0
  
  p1 <- exp(eta3)/(1+exp(eta3))
  p2 <- exp(eta2_y1*s2)/(1+exp(eta2_y1))
  p3 <- 1/(1+exp(eta3))
  p4 <- exp(eta2_y0*s2)/(1+exp(eta2_y0))
  l3 <- log(p1*p2 + p3*p4)
  
  lglik3 <- sum(l3)
  
  #loglikelihood
  loglik <- lglik1 + lglik2 + lglik3
  
  -loglik
  
}  


#####================================ WGLM ======================#############
mi<-10   #number of data sets from multiple imputation 
#Calculate weighting probability #ordinary logistic for Y|S,X in the validated data


lmod.svdata <- glm(y ~  s+crtdose1+crtdose2+crtdose3+hrtdose1 +hrtdose2 +hrtdose3+hypten+age, family=binomial(link="logit"), data=vdata)

coef.theta<- as.numeric(summary(lmod.svdata)$coefficients[, 1])
SE.theta<-as.numeric(summary(lmod.svdata)$coefficients[, 2])
var.theta<-vcov(lmod.svdata)
theta.M<-mvtnorm::rmvnorm(mi, mean = coef.theta, sigma =var.theta)


# Calculate the weight for nonvalidated sample 
#each individual has two weights == Calculate predicted P(Y=1|S_i,Z_i )and  P(Y=0|S_i,Z_i )==

cov.3 <- c("s" ,"crtdose1", "crtdose2", "crtdose3",  "hrtdose1","hrtdose2","hrtdose3", "hypten","age")


nvdat.y1 <- nvdata
nvdat.y1$y <- 1

mn.sy1 <- cbind(intercept=1, nvdat.y1[, cov.3])


### y=0
nvdat.y0 <- nvdata
nvdat.y0$y <- 0
mn.sy2 <- cbind(intercept=1, nvdat.y0[, cov.3])


coef.M<-matrix(rep(0,mi*6), nrow = mi, ncol = 6, byrow = TRUE)  #matrix of 10 coef.
SE.M<-matrix(rep(0,mi*6), nrow = mi, ncol = 6, byrow = TRUE)  #matrix of 10 SE


for (i in 1:mi){
  
mn.sy1<-as.matrix(mn.sy1)
eta1 <- mn.sy1 %*% theta.M[i,]

p.ygiven.sx <- exp(eta1)/(1+exp(eta1))


#combine vdata and nvdata with the weights as it is explained in the note for weighted GEE

modelvdata<-vdata
modelvdata$w<-1
modelnvdata1<-nvdata
modelnvdata2<-nvdata
modelnvdata1$y<-1
modelnvdata2$y<-0

modelnvdata1$w<-p.ygiven.sx
modelnvdata2$w<-1-p.ygiven.sx

modeldata<-rbind.fill(modelvdata,modelnvdata1,modelnvdata2)


#And now WGLM:#
modelWeightedGLM<-glm(y~hrtdose1 +hrtdose2 +hrtdose3+hypten+age ,  data=modeldata, na.action=na.omit,weights = w, family=binomial (link="logit") )
coef.M[i,]<- as.numeric(summary(modelWeightedGLM)$coefficients[, 1])
SE.M[i,]<-as.numeric(summary(modelWeightedGLM)$coefficients[, 2])
}


#======== Pooled estimate and standard error of M many imputation data set ======
coef<-colMeans(coef.M)


#==== To have a within and between variance I need variance, I need square of std

VW<-colMeans(SE.M^2)
VB<-apply(coef.M,2,var)
  
VTotal<-VW+(1+1/mi)*VB

SE<-sqrt(VTotal)

lb<-coef - 1.96 * SE
ub<-coef + 1.96 * SE

lb1<-round(lb,1)
ub1<-round(ub,1)


#######===================== Estimate the MLEs Using Pepe's method ==============###############
init1 <- rep(0,16)
ores1 <- optim(par=init1, fn=loglik.fnPepe, data=odata, method= 'BFGS', control = list(maxit=3000), hessian = TRUE) 
#ores1[["par"]]

hessian <- ores1$hessian
hessian.inv <- solve(hessian)
#hessian.inv <- ginv(hessian)
parameter.se <- sqrt(diag(hessian.inv))

lb<-ores1[["par"]][1:6] - 1.96 * parameter.se[1:6]
ub<-ores1[["par"]] [1:6]+ 1.96 * parameter.se[1:6]


lb2<-round(lb,1)
ub2<-round(ub,1)


################# logistic regression using self-reported outcome:#
modelCC<-glm(s~hrtdose1 +hrtdose2 +hrtdose3+hypten+age ,  data=tot.data, na.action=na.omit, family=binomial (link="logit") )
coefcc<- as.numeric(summary(modelCC)$coefficients[, 1])
SEcc<- as.numeric(summary(modelCC)$coefficients[, 2])


lb<-coefcc - 1.96 * SEcc
ub<-coefcc + 1.96 * SEcc


lb3<-round(lb,1)
ub3<-round(ub,1)


################# logistic regression using validated sub-sample:#
modelSJ<-glm(y~hrtdose1 +hrtdose2 +hrtdose3+hypten+age ,  data=vdata, na.action=na.omit, family=binomial (link="logit") )
#summary(modelSJLIFE)
coefsj<- as.numeric(summary(modelSJ)$coefficients[, 1])
SEsj<- as.numeric(summary(modelSJ)$coefficients[, 2])

lb<-coefsj - 1.96 * SEsj
ub<-coefsj + 1.96 * SEsj

lb4<-round(lb,1)
ub4<-round(ub,1)


#results<-c(ores1[["par"]], parameter.se ,coef,SE, coefcc,SEcc, coefsj,SEsj, V_Y.NV_N, V_N.NV_Y)

