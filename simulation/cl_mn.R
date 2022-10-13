# SIMULATION STUDY #
# Classical matching multinomial method #
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
  
  ################## CLASSICAL MATCHING MULTINOMIAL ##########################################
  
  tryCatch({
    # anatomic/semi-anatomic
    data_as<-data_sim[data_sim$treatment=="anatomic" | data_sim$treatment=="semi-anatomic",]
    data_as$treatment<-droplevels(data_as$treatment)
    match <- matchit(treatment ~ eta+creat+cirrosi+steato+
                       inv+MVI+HBV+HCV+satellitosi+sessoM, data=data_as,
                     method = "nearest", caliper=0.1)
    data_complete_match <- match.data(match)
    
    m51 <- coxph(Surv(T.event,Event)~ treatment+strata(subclass),data=data_complete_match)
    
    tab[j,41] <- summary(m51)$coef[2]
    tab[j,42] <- summary(m51)$coef[3]
    tab[j,43] <- summary(m51)$conf.int[3]
    tab[j,44] <- summary(m51)$conf.int[4]
    
    # anatomic/wedge
    data_aw<-data_sim[data_sim$treatment=="anatomic" | data_sim$treatment=="Wedge",]
    data_aw$treatment<-droplevels(data_aw$treatment)
    match <- matchit(treatment ~ eta+creat+cirrosi+steato+
                       inv+MVI+HBV+HCV+satellitosi+sessoM, data=data_aw,
                     method = "nearest", caliper=0.1)
    data_complete_match <- match.data(match)
    
    m52 <- coxph(Surv(T.event,Event)~ treatment+strata(subclass),data=data_complete_match)
    
    tab[j,45] <- summary(m52)$coef[2]
    tab[j,46] <- summary(m52)$coef[3]
    tab[j,47] <- summary(m52)$conf.int[3]
    tab[j,48] <- summary(m52)$conf.int[4]
    
    # Save tab
    write.csv(tab, file="sum.csv", row.names = FALSE)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
  print(j)
}