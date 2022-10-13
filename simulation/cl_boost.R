# SIMULATION STUDY #
# Classical matching boosting method #
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
  
  ########################### CLASSICAL MATCHING BOOSTING ####################################################
  
  tryCatch({
    data_as<-data_sim[data_sim$treatment=="anatomic" | data_sim$treatment=="semi-anatomic",]
    data_as$treatment<-droplevels(data_as$treatment)
    data_as$treatment<-ifelse(data_as$treatment=='semi-anatomic',1,0)
    match <- matchit(treatment ~ eta+creat+cirrosi+steato+
                       inv+MVI+HBV+HCV+satellitosi+sessoM, data=data_as,
                     distance="gbm", caliper = 0.1)
    data_complete_match <- match.data(match)
    
    m61 <- coxph(Surv(T.event,Event)~ treatment+strata(subclass),data=data_complete_match)
    
    tab[j,49] <- summary(m61)$coef[2]
    tab[j,50] <- summary(m61)$coef[3]
    tab[j,51] <- summary(m61)$conf.int[3]
    tab[j,52] <- summary(m61)$conf.int[4]
    
    # anatomic/wedge
    data_aw<-data_sim[data_sim$treatment=="anatomic" | data_sim$treatment=="Wedge",]
    data_aw$treatment<-droplevels(data_aw$treatment)
    data_aw$treatment<-ifelse(data_aw$treatment=="Wedge",1,0)
    match <- matchit(treatment ~ eta+creat+cirrosi+steato+
                       inv+MVI+HBV+HCV+satellitosi+sessoM, data=data_aw,
                     distance="gbm", caliper = 0.1)
    data_complete_match <- match.data(match)
    
    m62 <- coxph(Surv(T.event,Event)~ treatment+strata(subclass),data=data_complete_match)
    
    tab[j,53] <- summary(m62)$coef[2]
    tab[j,54] <- summary(m62)$coef[3]
    tab[j,55] <- summary(m62)$conf.int[3]
    tab[j,56] <- summary(m62)$conf.int[4]
    
    # Save tab
    write.csv(tab, file="sum.csv", row.names = FALSE)
    print(j)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
}