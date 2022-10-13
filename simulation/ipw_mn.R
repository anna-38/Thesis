# SIMULATION STUDY #
# Adjusted method #
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

N=1000

tab<-read.csv('sum.csv')

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
  
  ########################### IPW MULTINOMIAL ####################################################
  tryCatch({
    weights_mn <- weightit(treatment ~ eta+creat+cirrosi+steato+inv+
                             MVI+HBV+HCV+satellitosi+sessoM, data=data_sim, stabilize=TRUE, estimand="ATE")
    m1 <- coxph(Surv(T.event,Event)~ treatment,data=data_sim,weights = weights_mn$weights, robust = TRUE)
    
    tab[j,9] <- summary(m1)$coef[1,2]
    tab[j,10] <- summary(m1)$coef[1,3]
    tab[j,11] <- summary(m1)$conf.int[1,3]
    tab[j,12] <- summary(m1)$conf.int[1,4]
    
    tab[j,13] <- summary(m1)$coef[2,2]
    tab[j,14] <- summary(m1)$coef[2,3]
    tab[j,15] <- summary(m1)$conf.int[2,3]
    tab[j,16] <- summary(m1)$conf.int[2,4]
    # Save tab
    write.csv(tab, file="sum.csv", row.names = FALSE)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
  print(j)
}