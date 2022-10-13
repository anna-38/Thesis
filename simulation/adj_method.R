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

colu <- c('s_HR_adj','s_SE_adj','s_LOWER_adj','s_UPPER_adj','w_HR_adj','w_SE_adj','w_LOWER_adj','w_UPPER_adj',
          's_HR_ipw_mn','s_SE_ipw_mn','s_LOWER_ipw_mn','s_UPPER_ipw_mn','w_HR_ipw_mn','w_SE_ipw_mn','w_LOWER_ipw_mn','w_UPPER_ipw_mn',
          's_HR_ipw_boost','s_SE_ipw_boost','s_LOWER_ipw_boost','s_UPPER_ipw_boost','w_HR_ipw_boost','w_SE_ipw_boost','w_LOWER_ipw_boost','w_UPPER_ipw_boost',
          's_HR_mn_gps','s_SE_mn_gps','s_LOWER_mn_gps','s_UPPER_mn_gps','w_HR_mn_gps','w_SE_mn_gps','w_LOWER_mn_gps',
          'w_UPPER_mn_gps','s_HR_boost_gps','s_SE_boost_gps','s_LOWER_boost_gps','s_UPPER_boost_gps','w_HR_boost_gps',
          'w_SE_boost_gps','w_LOWER_boost_gps','w_UPPER_boost_gps','s_HR_mn_cl','s_SE_mn_cl','s_LOWER_mn_cl','s_UPPER_mn_cl','w_HR_mn_cl',
          'w_SE_mn_cl','w_LOWER_mn_cl','w_UPPER_mn_cl','s_HR_boost_cl','s_SE_boost_cl','s_LOWER_boost_cl','s_UPPER_boost_cl','w_HR_boost_cl','w_SE_boost_cl','w_LOWER_boost_cl','w_UPPER_boost_cl')

tab<-as.data.frame(matrix(nrow=1000,ncol=length(colu)))
colnames(tab)<-colu

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
  
  ########################### ADJUSTED MODEL ####################################################
  tryCatch({
    # eta,creat,sessoM,cirrosi,steato,HBV,HCV,satellitosi,inv,MVI
    mod_all<-coxph(Surv(T.event,Event)~ treatment+eta+creat+cirrosi+steato+
                     inv+MVI+HBV+HCV+satellitosi+sessoM,data=data_sim)
    
    tab[j,1] <- summary(mod_all)$coef[1,2]
    tab[j,2] <- summary(mod_all)$coef[1,3]
    tab[j,3] <- summary(mod_all)$conf.int[1,3]
    tab[j,4] <- summary(mod_all)$conf.int[1,4]
    
    tab[j,5] <- summary(mod_all)$coef[2,2]
    tab[j,6] <- summary(mod_all)$coef[2,3]
    tab[j,7] <- summary(mod_all)$conf.int[2,3]
    tab[j,8] <- summary(mod_all)$conf.int[2,4]
    # Save tab with results
    write.csv(tab, file="sum.csv", row.names = FALSE)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")})
  print(j)
}