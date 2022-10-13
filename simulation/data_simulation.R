## THESIS ##
## Monte Carlo method: Simulation data ##

rm(list=ls())

# import libraries
library(xlsx)
library(ggplot2)
library(cobalt)
library(prodlim)
library(survival)
library(WeightIt)
library(MatchIt)
library(mlogit)
library(gbm)
library(GPSCDF)
library(caret)
library(arm)

invisible(utils::memory.limit(64000))

# set the path
setwd("name_of_the_path")

# read dataset
data <-read.csv("final_data.csv")

# factor
var<-c('sesso','Cirrosi','Steatoepatite','child_pugh_grade','inv_macro','HBV','HCV','Potus.','PROCEDURE_cat','EdmondsonGrading','MVI','satellitosi','R')

data[,var]<-lapply(data[,var], function(x) factor(x))

set.seed(12347)

# For each method save:
# marginal estimated HR
# standard error (robust)
# confidence interval

# Methods:
# adjusted model
# IPW multinom
# IPW boosting
# GPS CDF matching multinom
# GPS CDF matching boosting
# classic matching multinom
# classic matching boosting

options(max.print=1000000)

all_data<-read.csv('sim.csv')

Nsim=1000
salti=0

for (j in 1:Nsim) {
  N=1000
  
  # generate the variables and the "alphas":
  
  # età (rnorm)
  # sessoM (rbinom)
  # cirrosi (rbinom)
  # Steatoepatite (rbinom)
  # HBV (rbinom)
  # HCV (rbinom)
  # creatininemia (rnorm)
  # satellitosi (rbinom)
  # inv_macro (rbinom)
  # MVI (rbinom)

  eta   <-  rnorm(N, mean(data$eta_surg), sd=sqrt(var(data$eta_surg)))
  alphaes <- 0.02 # semi vs anatomic
  alphaew <- 0.005 # wedge vs anatomic
  
  creat   <-  rnorm(N, mean(data$Creatininemia), sqrt(var(data$Creatininemia)))
  alphacs <- -0.04
  alphacw <- 0.05
  
  sessoM  <-  rbinom(N, 1, prob=summary(data$sesso)['M']/(summary(data$sesso)['M']+summary(data$sesso)['F']))
  alphass <- 0.05
  alphasw <- 0.01
  
  cirrosi  <-  rbinom(N, 1, prob=summary(data$Cirrosi)['1']/(summary(data$Cirrosi)['0']+summary(data$Cirrosi)['1']))
  alphacis <- 0.096
  alphaciw <- 0.005
  
  steato  <-  rbinom(N, 1, prob=summary(data$Steatoepatite)['1']/(summary(data$Steatoepatite)['0']+summary(data$Steatoepatite)['1']))
  alphasts <- 0.10
  alphastw <- -0.058
  
  HBV  <-  rbinom(N, 1, prob=summary(data$HBV)['1']/(summary(data$HBV)['0']+summary(data$HBV)['1']))
  alphahbs <- 0.079
  alphahbw <- -0.019
  
  HCV  <-  rbinom(N, 1, prob=summary(data$HCV)['1']/(summary(data$HCV)['0']+summary(data$HCV)['1']))
  alphahcs <- -0.01
  alphahcw <- -0.015
  
  satellitosi  <-  rbinom(N, 1, prob=summary(data$satellitosi)['1']/(summary(data$satellitosi)['0']+summary(data$satellitosi)['1']))
  alphasas <- 0.017
  alphasaw <- -0.025
  
  inv  <-  rbinom(N, 1, prob=summary(data$inv_macro)['1']/(summary(data$inv_macro)['0']+summary(data$inv_macro)['1']))
  alphais <- -1.28
  alphaiw <- -2.46
  
  MVI  <-  rbinom(N, 1, prob=summary(data$MVI)['1']/(summary(data$MVI)['0']+summary(data$MVI)['1']))
  alphams <- -1.26
  alphamw <- -0.082
  
  # ...
  
  alpha0s <- -1 # intercept semi vs anatomic (baseline hazard)
  alpha0w <- -0.9 # intercept wedged vs anatomic
  
  # probability of assigning the treatment anatomic, semi-anatomic or wedge
  ps <- invlogit(alpha0s+alphaes*eta+alphacs*creat+alphass*sessoM+alphacis*cirrosi+alphasts*steato+alphahbs*HBV+alphahcs*HCV+alphasas*satellitosi+alphais*inv+alphams*MVI)
  pw <- invlogit(alpha0w+alphaew*eta+alphacw*creat+alphasw*sessoM+alphaciw*cirrosi+alphastw*steato+alphahbw*HBV+alphahcw*HCV+alphasaw*satellitosi+alphaiw*inv+alphamw*MVI)
  
  treatment<-c()
  for( i in 1:length(ps)  ){
    if(ps[i]+pw[i]>=1){
      treatment.prob<-sample.int(n=3,size=1,prob=c(0,ps[i]/(ps[i]+pw[i]),pw[i]/(ps[i]+pw[i])))
    } else {
      treatment.prob<-sample.int(n=3,size=1,prob=c(1-ps[i]-pw[i],ps[i],pw[i]))
    }
    treatment<-c(treatment,ifelse(treatment.prob==1,"anatomic",ifelse(treatment.prob==2,"semi-anatomic","Wedged")))
  }
  
  t.cens <- runif(N, min=96, max=800) # in this way, no censoring
  U<-runif(N, min=0, max=1)
  k<- 0.0035 # Weibull scale
  p<- 1.5 # Weibull shape
  
  # beta inspired by multivariate Cox coefficients
  
  betae<- 0.01
  betac<- 0.49
  betas<- 0.11
  betaci<- 0.22
  betast<- -0.26
  betahb<- -0.02
  betahc<- -0.01
  betasa<- 0.34
  betai<- 0.48
  betam<- 0.52
  
  # Marginal hazard ratios
  HRs <-  0.2 # 0.2 # 1.5 # 1.5 # 0.8
  HRw <-  0.7 # 0.7 # 2 # 0.5 # 0.5
  
  beta0s <- log(HRs)
  beta0w <- log(HRw)
  iter = 0
  
  repeat {
    iter = iter+1
    # inverse of Survival Weibull: estimation of time of event
    t.event.a<- (-log(U))/(k*exp(betae*eta+betac*creat+betas*sessoM+betaci*cirrosi+betast*steato+betahb*HBV+betahc*HCV+betasa*satellitosi+betai*inv+betam*MVI))^(1/p)  
    t.event.s <- (-log(U))/(k*exp(beta0s+betae*eta+betac*creat+betas*sessoM+betaci*cirrosi+betast*steato+betahb*HBV+betahc*HCV+betasa*satellitosi+betai*inv+betam*MVI))^(1/p)  
    mod <- coxph(Surv(c(t.event.a,t.event.s),rep(1,2*N))~c(rep(0,N),rep(1,N)))
    if (round(summary(mod)$coef[1,1],3) == round(log(HRs),3) | iter==1000) break
    if (summary(mod)$coef[1,1] < log(HRs)) {beta0s <- beta0s + 0.001}
    if (summary(mod)$coef[1,1] > log(HRs)) {beta0s <- beta0s - 0.001}
    
  }
  
  if (iter ==1000) {salti = salti+1 # drop the dataset if iterations>1000
  next
  }
  
  repeat {
    iter = iter+1
    t.event.a<- (-log(U))/(k*exp(betae*eta+betac*creat+betas*sessoM+betaci*cirrosi+betast*steato+betahb*HBV+betahc*HCV+betasa*satellitosi+betai*inv+betam*MVI))^(1/p)  
    t.event.w <- (-log(U))/(k*exp(beta0w+betae*eta+betac*creat+betas*sessoM+betaci*cirrosi+betast*steato+betahb*HBV+betahc*HCV+betasa*satellitosi+betai*inv+betam*MVI))^(1/p)  
    mod <- coxph(Surv(c(t.event.a,t.event.w),rep(1,2*N))~c(rep(0,N),rep(1,N)))
    if (round(summary(mod)$coef[1,1],3) == round(log(HRw),3) | iter==1000) break
    if (summary(mod)$coef[1,1] < log(HRw)) {beta0w <- beta0w + 0.001}
    if (summary(mod)$coef[1,1] > log(HRw)) {beta0w <- beta0w - 0.001}
    
  }
  
  if (iter ==1000) {
    salti = salti+1
    next
  }
  
  
  
  # variable of the time of event
  # and binary variable that indicates if the event happened or not
  t.event <- ifelse(treatment=="anatomic", t.event.a, ifelse(treatment=="semi-anatomic",t.event.s,t.event.w))
  T.event<-pmin(t.event,t.cens)
  Event<-ifelse(t.event<t.cens,1,0)
  
  # CREATE SIMULATED DATASET
  data_sim <- data.frame(T.event,Event,treatment, eta,creat,sessoM,cirrosi,steato,HBV,HCV,satellitosi,inv,MVI)
  var<-c('treatment','sessoM','cirrosi','steato','HBV','HCV','satellitosi','inv','MVI')
  data_sim[,var]<-lapply(data_sim[,var], function(x) factor(x))
  
  if(j==1){
    all_data<-data_sim
  }else{
    all_data<-rbind(all_data,data_sim)
  }
  write.csv(all_data, file='sim.csv', row.names = FALSE)
  print(j)
}

# sim_A HRs=0.8, HRw=0.5
# sim_B HRs=1.5, HRw=0.5
# sim_C HRs=1.5, HRw=2
# sim_D HRs=0.2, HRw=0.7
