# SIMULATION STUDY #
# GPS-CDF multinomial method #
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

# read dataset that contains the simulations
data <-read.csv("sim.csv")

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

set.seed(12347)
options(max.print=1000000)

tab<-read.csv('sum.csv')

N=1000
for (j in 1:N) {
  # SIMULATED DATASET
  
  data_sim<-data[c((1+(j-1)*1000):(j*1000)),]
  var<-c('treatment','sessoM','cirrosi','steato','HBV','HCV','satellitosi','inv','MVI')
  data_sim[,var]<-lapply(data_sim[,var], function(x) factor(x))
  #################################################
  # APPLICATION OF THE METHODS TO THE SIMULATED DATASETS
  #################################################
  
  # Methods:
  # adjusted model
  # IPW multinom
  # IPW boosting
  # GPS CDF matching multinom
  # GPS CDF matching boosting
  # classic matching multinom
  # classic matching boosting
  
  ########################### GPS-CDF MULTINOMIAL ####################################################
  tryCatch({
    adj_mod<-nnet::multinom(as.factor(treatment) ~ eta+creat+cirrosi+steato+
                              inv+MVI+HBV+HCV+satellitosi+sessoM, data=data_sim, verbose = FALSE)
    probab<- round(predict(adj_mod, newdata=data_sim, type="prob"),digits=8)
    gps<-cbind(probab[,1],probab[,2],1-probab[,1]-probab[,2])
    proc<-as.numeric(as.factor(data_sim$treatment))
    fit<-GPSCDF(pscores=gps, data=data_sim, trt=proc, greedy=TRUE, ordinal=TRUE)
    
    overall_sets<-table(fit$grddata$treatment, fit$grdmatch)
    ind_as<-as.numeric(colnames(overall_sets[,overall_sets[3,]==0]))
    ind_aw<-as.numeric(colnames(overall_sets[,overall_sets[2,]==0]))
    
    # anatomic/semi-anatomic
    data_tmp<-fit$grddata[fit$grddata$grdmatch %in% ind_as,]
    data_tmp$treatment<-droplevels(data_tmp$treatment)
    m31 <- coxph(Surv(T.event,Event)~ treatment+strata(grdmatch),data=data_tmp)
    
    tab[j,25] <- summary(m31)$coef[2]
    tab[j,26] <- summary(m31)$coef[3]
    tab[j,27] <- summary(m31)$conf.int[3]
    tab[j,28] <- summary(m31)$conf.int[4]
    
    # semi-anatomic/wedge
    data_tmp<-fit$grddata[fit$grddata$grdmatch %in% ind_aw,]
    data_tmp$treatment<-droplevels(data_tmp$treatment)
    m32 <- coxph(Surv(T.event,Event)~ treatment+strata(grdmatch),data=data_tmp)
    
    tab[j,29] <- summary(m32)$coef[2]
    tab[j,30] <- summary(m32)$coef[3]
    tab[j,31] <- summary(m32)$conf.int[3]
    tab[j,32] <- summary(m32)$conf.int[4]
    
    # Save tab
    write.csv(tab, file="sum.csv", row.names = FALSE)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
  print(j)
}