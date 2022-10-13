# SIMULATION STUDY #
# IPW boosting method #
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
  # DATASET SIMULATO
  
  data_sim<-data[c((1+(j-1)*N):(j*N)),]
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
  
  ########################### IPW BOOSTING ####################################################
  tryCatch({
    weights_gbm <- weightit(treatment ~ eta+creat+cirrosi+steato+
                              inv+MVI+HBV+HCV+satellitosi+sessoM, data=data_sim,
                            method = "gbm", stabilize=TRUE, stop.method = 'es.mean')
    m2 <- coxph(Surv(T.event,Event)~ treatment,data=data_sim,weights = weights_gbm$weights, robust = TRUE)
    
    tab[j,17] <- summary(m2)$coef[1,2]
    tab[j,18] <- summary(m2)$coef[1,3]
    tab[j,19] <- summary(m2)$conf.int[1,3]
    tab[j,20] <- summary(m2)$conf.int[1,4]
    
    tab[j,21] <- summary(m2)$coef[2,2]
    tab[j,22] <- summary(m2)$coef[2,3]
    tab[j,23] <- summary(m2)$conf.int[2,3]
    tab[j,24] <- summary(m2)$conf.int[2,4]
    # Save tab
    write.csv(tab, file="sum.csv", row.names = FALSE)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
  print(j)
}