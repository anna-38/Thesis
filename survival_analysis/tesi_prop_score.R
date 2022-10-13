## THESIS ##
## Propensity score method ##
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
library(survminer)
library(caret)

# set the path
setwd("name_of_the_path")

# read dataset
data <-read.csv("final_data.csv")

# factor
var<-c('sesso','Cirrosi','Steatoepatite','child_pugh_grade','inv_macro','HBV','HCV','Potus.','PROCEDURE_cat','EdmondsonGrading','MVI','satellitosi','R')

data$PROCEDURE_cat<-as.factor(data$PROCEDURE_cat)
data$sesso<-as.factor(data$sesso)
data[,var]<-lapply(data[,var], function(x) factor(x))

# tables to save the data
colu <- c('s_HR','s_SE','s_LOWER','s_UPPER','s_pvalue','w_HR','w_SE','w_LOWER','w_UPPER','w_pvalue')
tabOS<-as.data.frame(matrix(nrow=6,ncol=length(colu)))
colnames(tabOS)<-colu
row.names(tabOS)<-c('unadjusted','adjusted','mult_weights','boost_weights','mult_gps','mult_classic')

tabRFS<-as.data.frame(matrix(nrow=6,ncol=length(colu)))
colnames(tabRFS)<-colu
row.names(tabRFS)<-c('unadjusted','adjusted','mult_weights','boost_weights','mult_gps','mult_classic')

################# COX MODELS ##################################################################################################################################################################################################

## UNADJUSTED MODEL ##
##  Univariate Cox  ##

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="bottomright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
title('OS unadjusted model')
summary(sFit,times = c(12,24,36,48,60,120))

# OS

cox<- coxph(Surv(timeOSmesi,eventomorte)~PROCEDURE_cat,data)

summary(cox)
tabOS[1,1] <- summary(cox)$coef[1,2]
tabOS[1,2] <- summary(cox)$coef[1,3]
tabOS[1,3] <- summary(cox)$conf.int[1,3]
tabOS[1,4] <- summary(cox)$conf.int[1,4]
tabOS[1,5] <- summary(cox)$coef[1,5]

tabOS[1,6] <- summary(cox)$coef[2,2]
tabOS[1,7] <- summary(cox)$coef[2,3]
tabOS[1,8] <- summary(cox)$conf.int[2,3]
tabOS[1,9] <- summary(cox)$conf.int[2,4]
tabOS[1,10] <- summary(cox)$coef[2,5]

# RFS

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.title='Resection',atrisk.cex=0.8,legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
title('RFS unadjusted model')
summary(sFit,times = c(12,24,36,48,60,120))

cox<- coxph(Surv(timeRFSmesi,RFS)~PROCEDURE_cat, data)

summary(cox)
tabRFS[1,1] <- summary(cox)$coef[1,2]
tabRFS[1,2] <- summary(cox)$coef[1,3]
tabRFS[1,3] <- summary(cox)$conf.int[1,3]
tabRFS[1,4] <- summary(cox)$conf.int[1,4]
tabRFS[1,5] <- summary(cox)$coef[1,5]

tabRFS[1,6] <- summary(cox)$coef[2,2]
tabRFS[1,7] <- summary(cox)$coef[2,3]
tabRFS[1,8] <- summary(cox)$conf.int[2,3]
tabRFS[1,9] <- summary(cox)$conf.int[2,4]
tabRFS[1,10] <- summary(cox)$coef[2,5]

## ADJUSTED MODEL ##
## Cox model with all variables ##

# OS

cox_all<- coxph(Surv(timeOSmesi,eventomorte)~eta_surg+bili_tot+Albuminemia+Creatininemia+
                  INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                  MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso+PROCEDURE_cat,data)

summary(cox_all)
tabOS[2,1] <- summary(cox_all)$coef[21,2]
tabOS[2,2] <- summary(cox_all)$coef[21,3]
tabOS[2,3] <- summary(cox_all)$conf.int[21,3]
tabOS[2,4] <- summary(cox_all)$conf.int[21,4]
tabOS[2,5] <- summary(cox_all)$coef[21,5]

tabOS[2,6] <- summary(cox_all)$coef[22,2]
tabOS[2,7] <- summary(cox_all)$coef[22,3]
tabOS[2,8] <- summary(cox_all)$conf.int[22,3]
tabOS[2,9] <- summary(cox_all)$conf.int[22,4]
tabOS[2,10] <- summary(cox_all)$coef[22,5]

# RFS

cox_all<- coxph(Surv(timeRFSmesi,RFS)~eta_surg+bili_tot+Albuminemia+Creatininemia+
                  INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                  MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso+PROCEDURE_cat, data)

summary(cox_all)
tabRFS[2,1] <- summary(cox_all)$coef[21,2]
tabRFS[2,2] <- summary(cox_all)$coef[21,3]
tabRFS[2,3] <- summary(cox_all)$conf.int[21,3]
tabRFS[2,4] <- summary(cox_all)$conf.int[21,4]
tabRFS[2,5] <- summary(cox_all)$coef[21,5]

tabRFS[2,6] <- summary(cox_all)$coef[22,2]
tabRFS[2,7] <- summary(cox_all)$coef[22,3]
tabRFS[2,8] <- summary(cox_all)$conf.int[22,3]
tabRFS[2,9] <- summary(cox_all)$conf.int[22,4]
tabRFS[2,10] <- summary(cox_all)$coef[22,5]

#####################################################################################
########################## PROPENSITY SCORE METHODS ##################################
### IPW multinomial ###
### IPW boosting ###
### IPW boosting ###
### GPS-CDF multinomial ###
### GPS-CDF boosting ###
### Classical matching multinomial ###
### Classical matching boosting ###


########################## Propensity score weights IPW ###############################
## For ATE ##

################################# MULTINOMIAL ##########################################
## SELECTED VARIABLES: ##
# eta_surg+bili_tot+Albuminemia+Creatininemia+
# INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
# MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso

weights_mn <- weightit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                         INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                         MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                       data=data, stabilize=TRUE) # method = "ps" default: corresponds to multinomial
weights_mn # weights computed
summary(weights_mn)

# Balance plots:

# - overall
par(mfrow=c(1,1),mar=c(4,4,1,1))
boxplot(weights_mn$weights~PROCEDURE_cat,data=data, col=c('orange','violet','lightblue'), xaxt='n', ylab="", xlab='resection')
axis(1, at=c(1,2,3), labels = FALSE)
text(x =  seq_along(c('anatomic','semi-anatomic','wedge')), y = -0.5, srt = 20, adj = 1,
     labels = c('anatomic','semi-anatomic','wedge'), xpd = TRUE)
title('IPW multinomial weights')

# - anatomic/semi-anatomic
weights_as<-weights_mn
weights_as$weights<-weights_mn$weights[weights_mn$treat=="anatomic" | weights_mn$treat=="semi-anatomic"]
weights_as$treat<-weights_mn$treat[weights_mn$treat=="anatomic" | weights_mn$treat=="semi-anatomic"]
weights_as$covs<-weights_mn$covs[weights_mn$treat=="anatomic" | weights_mn$treat=="semi-anatomic",]
weights_as$s.weights<-weights_mn$s.weights[weights_mn$treat=="anatomic" | weights_mn$treat=="semi-anatomic"]

love.plot(weights_as,threshold = .1,drop.distance = T,
          stars="raw",stats="m",continuous="std",binary="std",
          var.names = data.frame(names(weights_mn$covs)), title = "Covariate Balance anatomic vs semi-anatomic")
bal.tab(weights_as, stats = "m", thresholds = .1, binary="std")

# - semi-anatomic/wedge
weights_sw<-weights_mn
weights_sw$weights<-weights_mn$weights[weights_mn$treat=="semi-anatomic" | weights_mn$treat=="Wedge"]
weights_sw$treat<-weights_mn$treat[weights_mn$treat=="semi-anatomic" | weights_mn$treat=="Wedge"]
weights_sw$covs<-weights_mn$covs[weights_mn$treat=="semi-anatomic" | weights_mn$treat=="Wedge",]
weights_sw$s.weights<-weights_mn$s.weights[weights_mn$treat=="semi-anatomic" | weights_mn$treat=="Wedge"]

love.plot(weights_sw,threshold = .1,drop.distance = T,
          stars="raw",stats="m",continuous="std",binary="std",
          var.names = data.frame(names(weights_mn$covs)))
bal.tab(weights_sw, stats = "m", thresholds = .1, binary="std")

# - anatomic/wedge
weights_aw<-weights_mn
weights_aw$weights<-weights_mn$weights[weights_mn$treat=="anatomic" | weights_mn$treat=="Wedge"]
weights_aw$treat<-weights_mn$treat[weights_mn$treat=="anatomic" | weights_mn$treat=="Wedge"]
weights_aw$covs<-weights_mn$covs[weights_mn$treat=="anatomic" | weights_mn$treat=="Wedge",]
weights_aw$s.weights<-weights_mn$s.weights[weights_mn$treat=="anatomic" | weights_mn$treat=="Wedge"]

love.plot(weights_aw,threshold = .1,drop.distance = T,
          stars="raw",stats="m",continuous="std",binary="std",
          var.names = data.frame(names(weights_mn$covs)), title = "Covariate Balance anatomic vs wedge")
bal.tab(weights_aw, stats = "m", thresholds = .1, binary="std")

## MODELS

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data, caseweights = weights_mn$weights)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="bottomright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS ipw multinomial')

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat,data=data,weights = weights_mn$weights, robust = TRUE)
summary(m)

tabOS[3,1] <- summary(m)$coef[1,2]
tabOS[3,2] <- summary(m)$coef[1,3]
tabOS[3,3] <- summary(m)$conf.int[1,3]
tabOS[3,4] <- summary(m)$conf.int[1,4]
tabOS[3,5] <- summary(m)$coef[1,6]

tabOS[3,6] <- summary(m)$coef[2,2]
tabOS[3,7] <- summary(m)$coef[2,3]
tabOS[3,8] <- summary(m)$conf.int[2,3]
tabOS[3,9] <- summary(m)$conf.int[2,4]
tabOS[3,10] <- summary(m)$coef[2,6]

# RFS

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data, caseweights = weights_mn$weights)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('RFS ipw multinomial')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat,data=data,weights = weights_mn$weights, robust = TRUE)
summary(m)

tabRFS[3,1] <- summary(m)$coef[1,2]
tabRFS[3,2] <- summary(m)$coef[1,3]
tabRFS[3,3] <- summary(m)$conf.int[1,3]
tabRFS[3,4] <- summary(m)$conf.int[1,4]
tabRFS[3,5] <- summary(m)$coef[1,6]

tabRFS[3,6] <- summary(m)$coef[2,2]
tabRFS[3,7] <- summary(m)$coef[2,3]
tabRFS[3,8] <- summary(m)$conf.int[2,3]
tabRFS[3,9] <- summary(m)$conf.int[2,4]
tabRFS[3,10] <- summary(m)$coef[2,6]

################################# BOOSTING ##########################################

weights_gbm <- weightit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                          INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                          MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                        data=data, method = "gbm", stabilize=TRUE, stop.method = 'es.mean') # method = "gbm" boosting, stop method: the suggested one (es=standardized effect size)
weights_gbm # weights computed
summary(weights_gbm)

# Balance plots:

# - overall
par(mfrow=c(1,1),mar=c(4,4,1,1))
boxplot(weights_gbm$weights~PROCEDURE_cat,data=data, col=c('orange','violet','lightblue'), xaxt='n', ylab="", xlab='resection')
axis(1, at=c(1,2,3), labels = FALSE)
text(x =  seq_along(c('anatomic','semi-anatomic','wedge')), y=0.1, srt = 20, adj = 1, #y = -0.5, 
     labels = c('anatomic','semi-anatomic','wedge'), xpd = TRUE)
title('IPW gbm weights')

# - anatomic/semi-anatomic
weights_as<-weights_gbm
weights_as$weights<-weights_gbm$weights[weights_gbm$treat=="anatomic" | weights_gbm$treat=="semi-anatomic"]
weights_as$treat<-weights_gbm$treat[weights_gbm$treat=="anatomic" | weights_gbm$treat=="semi-anatomic"]
weights_as$covs<-weights_gbm$covs[weights_gbm$treat=="anatomic" | weights_gbm$treat=="semi-anatomic",]
weights_as$s.weights<-weights_gbm$s.weights[weights_gbm$treat=="anatomic" | weights_gbm$treat=="semi-anatomic"]

love.plot(weights_as,threshold = .1,drop.distance = T,
          stars="raw",stats="m",continuous="std",binary="std",
          var.names = data.frame(names(weights_gbm$covs)), title = "Covariate Balance anatomic vs semi-anatomic")
bal.tab(weights_as, stats = "m", thresholds = .1, binary="std")

# - semi-anatomic/wedge
weights_sw<-weights_gbm
weights_sw$weights<-weights_gbm$weights[weights_gbm$treat=="semi-anatomic" | weights_gbm$treat=="Wedge"]
weights_sw$treat<-weights_gbm$treat[weights_gbm$treat=="semi-anatomic" | weights_gbm$treat=="Wedge"]
weights_sw$covs<-weights_gbm$covs[weights_gbm$treat=="semi-anatomic" | weights_gbm$treat=="Wedge",]
weights_sw$s.weights<-weights_gbm$s.weights[weights_gbm$treat=="semi-anatomic" | weights_gbm$treat=="Wedge"]

love.plot(weights_sw,threshold = .1,drop.distance = T,
          stars="raw",stats="m",continuous="std",binary="std",
          var.names = data.frame(names(weights_gbm$covs)))
bal.tab(weights_sw, stats = "m", thresholds = .1, binary="std")

# - anatomic/wedge
weights_aw<-weights_gbm
weights_aw$weights<-weights_gbm$weights[weights_gbm$treat=="anatomic" | weights_gbm$treat=="Wedge"]
weights_aw$treat<-weights_gbm$treat[weights_gbm$treat=="anatomic" | weights_gbm$treat=="Wedge"]
weights_aw$covs<-weights_gbm$covs[weights_gbm$treat=="anatomic" | weights_gbm$treat=="Wedge",]
weights_aw$s.weights<-weights_gbm$s.weights[weights_gbm$treat=="anatomic" | weights_gbm$treat=="Wedge"]

love.plot(weights_aw,threshold = .1,drop.distance = T,
          stars="raw",stats="m",continuous="std",binary="std",
          var.names = data.frame(names(weights_mn$covs)), title = "Covariate Balance anatomic vs wedge") # anche qui inv_macro_sovraep molto sbilanciato
bal.tab(weights_aw, stats = "m", thresholds = .1, binary="std")

## MODELS

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data, caseweights = weights_gbm$weights)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="bottomright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS ipw gbm')

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat,data=data, weights = weights_gbm$weights, robust = TRUE)
summary(m)

tabOS[4,1] <- summary(m)$coef[1,2]
tabOS[4,2] <- summary(m)$coef[1,3]
tabOS[4,3] <- summary(m)$conf.int[1,3]
tabOS[4,4] <- summary(m)$conf.int[1,4]
tabOS[4,5] <- summary(m)$coef[1,6]

tabOS[4,6] <- summary(m)$coef[2,2]
tabOS[4,7] <- summary(m)$coef[2,3]
tabOS[4,8] <- summary(m)$conf.int[2,3]
tabOS[4,9] <- summary(m)$conf.int[2,4]
tabOS[4,10] <- summary(m)$coef[2,6]

# RFS

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data, caseweights = weights_gbm$weights)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Disease Free survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('RFS ipw gbm')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat,data=data,weights = weights_gbm$weights, robust = TRUE)
summary(m)

tabRFS[4,1] <- summary(m)$coef[1,2]
tabRFS[4,2] <- summary(m)$coef[1,3]
tabRFS[4,3] <- summary(m)$conf.int[1,3]
tabRFS[4,4] <- summary(m)$conf.int[1,4]
tabRFS[4,5] <- summary(m)$coef[1,6]

tabRFS[4,6] <- summary(m)$coef[2,2]
tabRFS[4,7] <- summary(m)$coef[2,3]
tabRFS[4,8] <- summary(m)$conf.int[2,3]
tabRFS[4,9] <- summary(m)$conf.int[2,4]
tabRFS[4,10] <- summary(m)$coef[2,6]

##################################### GPSCDF ##############################################

# Create the generalized propensity score (GPS) vector using any parametric or
# non-parametric model

## For ATE ##

############################### MULTINOMIAL ###################################

glm<- nnet::multinom(as.factor(PROCEDURE_cat) ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                       INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                       MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, data=data)
probab<- round(predict(glm, newdata=data, type="probs"),digits=8)
gps<-cbind(probab[,1],probab[,2],1-probab[,1]-probab[,2])
# Create scalar balancing power parameter
proc<-as.numeric(as.factor(data$PROCEDURE_cat))
fit<-GPSCDF(pscores=gps, data=data, trt=proc, optimal=TRUE, ordinal=TRUE)

# This method created a matched data, where one patient is associated to ONLY ONE patient (treated with another treatment).
# That means that there aren't sets of three patients treated with all the three treatments, but only sets with patients treated with 2 of the 3 procedures.
# So, now let's compare the graphics considering the three sets (anatomic/semi-anatomic; anatomic/wedge; wedge/semi-anatomic)

overall_sets<-table(fit$data$PROCEDURE_cat, fit$optmatch)
ind_as<-as.numeric(colnames(overall_sets[,overall_sets[3,]==0]))
ind_aw<-as.numeric(colnames(overall_sets[,overall_sets[2,]==0]))
ind_sw<-as.numeric(colnames(overall_sets[,overall_sets[1,]==0]))

# anatomic/semi-anatomic
data_tmp<-fit$data[fit$data$optmatch %in% ind_as,]
data_tmp$PROCEDURE_cat<-droplevels(data_tmp$PROCEDURE_cat)

covs<-data_tmp[,c('eta_surg','bili_tot','Albuminemia','Creatininemia',
        'INR','Piastrine','Cirrosi','child_pugh_grade','Steatoepatite','inv_macro',
        'MVI','EdmondsonGrading','HBV','HCV','Potus.','satellitosi','R','sesso')]

# covs_un<-data_as[,c('eta_surg','bili_tot','Albuminemia','Creatininemia',
#                   'INR','Piastrine','Cirrosi','child_pugh_grade','Steatoepatite','inv_macro',
#                   'MVI','EdmondsonGrading','HBV','HCV','Potus.','satellitosi','R','sesso')]

love.plot(bal.tab(PROCEDURE_cat ~ covs,data=data_tmp), threshold = .1,drop.distance = T,
        stars="raw",stats="m",abs=TRUE,continuous="std",binary="std", sample.names = 'Matched', colors = 'cyan3', limits=c(0,0.4), title = "Covariate Balance anatomic vs semi-anatomic")

# love.plot(bal.tab(PROCEDURE_cat ~ covs_un,data=data_as), threshold = .1,drop.distance = T,
#           stars="raw",stats="m",abs=TRUE,continuous="std",binary="std", title = "Covariate Balance anatomic vs semi-anatomic")

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS gps-cdf multinomial')

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(optmatch),data=data_tmp)
summary(m)

tabOS[5,6] <- summary(m)$coef[2]
tabOS[5,7] <- summary(m)$coef[3]
tabOS[5,8] <- summary(m)$conf.int[3]
tabOS[5,9] <- summary(m)$conf.int[4]
tabOS[5,10] <- summary(m)$coef[5]

# RFS

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('RFS gps-cdf multinomial')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(optmatch),data=data_tmp)
summary(m)

tabRFS[5,6] <- summary(m)$coef[2]
tabRFS[5,7] <- summary(m)$coef[3]
tabRFS[5,8] <- summary(m)$conf.int[3]
tabRFS[5,9] <- summary(m)$conf.int[4]
tabRFS[5,10] <- summary(m)$coef[5]

# semi-anatomic/wedge
data_tmp<-fit$grddata[fit$grddata$grdmatch %in% ind_sw,]
data_tmp$PROCEDURE_cat<-droplevels(data_tmp$PROCEDURE_cat)
sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS gps-cdf multinomial')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="bottomright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

# anatomic/wedge
data_tmp<-fit$data[fit$data$optmatch %in% ind_aw,]
data_tmp$PROCEDURE_cat<-droplevels(data_tmp$PROCEDURE_cat)

covs<-data_tmp[,c('eta_surg','bili_tot','Albuminemia','Creatininemia',
                  'INR','Piastrine','Cirrosi','child_pugh_grade','Steatoepatite','inv_macro',
                  'MVI','EdmondsonGrading','HBV','HCV','Potus.','satellitosi','R','sesso')]

# covs_un<-data_aw[,c('eta_surg','bili_tot','Albuminemia','Creatininemia',
#                     'INR','Piastrine','Cirrosi','child_pugh_grade','Steatoepatite','inv_macro',
#                     'MVI','EdmondsonGrading','HBV','HCV','Potus.','satellitosi','R','sesso')]
# love.plot(bal.tab(PROCEDURE_cat ~ covs,data=data_tmp), threshold = .1,drop.distance = T,
#           stars="raw",stats="m",abs=TRUE,continuous="std",binary="std", sample.names = 'Matched', colors = 'cyan3', limits=c(0,0.72), title = "Covariate Balance anatomic vs wedge")

love.plot(bal.tab(PROCEDURE_cat ~ covs_un,data=data_aw), threshold = .1,drop.distance = T,
          stars="raw",stats="m",abs=TRUE,continuous="std",binary="std", limits=c(0,0.72), title = "Covariate Balance anatomic vs semi-anatomic")

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS gps-cdf multinomial')

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(optmatch),data=data_tmp)
summary(m)

tabOS[5,6] <- summary(m)$coef[2]
tabOS[5,7] <- summary(m)$coef[3]
tabOS[5,8] <- summary(m)$conf.int[3]
tabOS[5,9] <- summary(m)$conf.int[4]
tabOS[5,10] <- summary(m)$coef[5]

# RFS
sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('RFS gps-cdf multinomial')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(optmatch),data=data_tmp)
summary(m)

tabRFS[5,6] <- summary(m)$coef[2]
tabRFS[5,7] <- summary(m)$coef[3]
tabRFS[5,8] <- summary(m)$conf.int[3]
tabRFS[5,9] <- summary(m)$conf.int[4]
tabRFS[5,10] <- summary(m)$coef[5]

################### CHECK ON ORDER ##############################################################
######### It's just a security test, not obligatory to run ######################################

# change order: should stay the same?
# No, because there is an order, but, if I reverse the order, the result is the same
proc2<-proc
proc2[proc==1]<-3
proc2[proc==3]<-1
fit2<-GPSCDF(pscores=gps, data=data, trt=proc, greedy=TRUE, ordinal=TRUE)

overall_sets<-table(fit2$grddata$PROCEDURE_cat, fit2$grdmatch)
ind_as<-as.numeric(colnames(overall_sets[,overall_sets[3,]==0]))
ind_aw<-as.numeric(colnames(overall_sets[,overall_sets[2,]==0]))
ind_sw<-as.numeric(colnames(overall_sets[,overall_sets[1,]==0]))

# anatomic/semi-anatomic
data_tmp<-fit2$grddata[fit2$grddata$grdmatch %in% ind_as,]
sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Disease Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

# The results are identical! As we expected

#############################################################################################

##################################### BOOSTING ##############################################
set.seed(123)
boost<- train(as.factor(PROCEDURE_cat) ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, data=data, method = "gbm",
              verbose = FALSE)
probab<- round(predict(boost, newdata=data, type="prob"),digits=8)
gps<-cbind(probab[,1],probab[,2],1-probab[,1]-probab[,2])
# Create scalar balancing power parameter
proc<-as.numeric(as.factor(data$PROCEDURE_cat))
fit<-GPSCDF(pscores=gps, data=data, trt=proc, greedy=TRUE, ordinal=TRUE)

overall_sets<-table(fit$grddata$PROCEDURE_cat, fit$grdmatch)
ind_as<-as.numeric(colnames(overall_sets[,overall_sets[3,]==0]))
ind_aw<-as.numeric(colnames(overall_sets[,overall_sets[2,]==0]))
ind_sw<-as.numeric(colnames(overall_sets[,overall_sets[1,]==0]))

# anatomic/semi-anatomic
data_tmp<-fit$grddata[fit$grddata$grdmatch %in% ind_as,]
data_tmp$PROCEDURE_cat<-droplevels(data_tmp$PROCEDURE_cat)
sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Disease Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

# semi-anatomic/wedge
data_tmp<-fit$grddata[fit$grddata$grdmatch %in% ind_sw,]
data_tmp$PROCEDURE_cat<-droplevels(data_tmp$PROCEDURE_cat)
sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Disease Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

# anatomic/wedge
data_tmp<-fit$grddata[fit$grddata$grdmatch %in% ind_aw,]
data_tmp$PROCEDURE_cat<-droplevels(data_tmp$PROCEDURE_cat)
sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Disease Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat, data=data_tmp)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(grdmatch),data=data_tmp)
summary(m)

##################################### CLASSICAL MODEL ##############################################
############### Consider models with binary outcomes and compare them ##############################
############### In particular, we use MATCHING PS METHOD ###########################################

## For ATT ##

##################################### MULTINOMIAL ##############################################

# anatomic/semi-anatomic
data_as<-data[data$PROCEDURE_cat=="anatomic" | data$PROCEDURE_cat=="semi-anatomic",]
data_as$PROCEDURE_cat<-droplevels(data_as$PROCEDURE_cat)
match <- matchit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                   INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                   MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                 data=data_as,
                 method = "nearest", caliper=0.3) # default distance="glm" generalized linear model
data_complete_match <- match.data(match)

# Balance plots:
love.plot(match,threshold = .1,drop.distance = T,sample.names=c("Unmatched", "Matched"),
          stars="raw",stats="m",continuous="std",binary="std", title='Covariate Balance Anatomic vs Semi-anatomic')

bal.tab(match, stats = "m", thresholds = .1, binary="std")

## MODELS

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="bottomright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS matching')

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

tabOS[6,1] <- summary(m)$coef[2]
tabOS[6,2] <- summary(m)$coef[3]
tabOS[6,3] <- summary(m)$conf.int[3]
tabOS[6,4] <- summary(m)$conf.int[4]
tabOS[6,5] <- summary(m)$coef[5]

#RFS

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('RFS matching')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

tabRFS[6,1] <- summary(m)$coef[2]
tabRFS[6,2] <- summary(m)$coef[3]
tabRFS[6,3] <- summary(m)$conf.int[3]
tabRFS[6,4] <- summary(m)$conf.int[4]
tabRFS[6,5] <- summary(m)$coef[5]

# anatomic/wedge
data_aw<-data[data$PROCEDURE_cat=="anatomic" | data$PROCEDURE_cat=="Wedge",]
data_aw$PROCEDURE_cat<-droplevels(data_aw$PROCEDURE_cat)
match <- matchit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                   INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                   MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                 data=data_aw,
                 method = "nearest", caliper=0.3)
data_complete_match <- match.data(match)

# Balance plots:

love.plot(match,threshold = .1,drop.distance = T,sample.names=c("Unmatched", "Matched"),
          stars="raw",stats="m",continuous="std",binary="std", title='Covariate Balance Anatomic vs Wedged')

bal.tab(match, stats = "m", thresholds = .1, binary="std")

## MODELS

# OS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="bottomright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('OS matching')

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

tabOS[6,6] <- summary(m)$coef[2]
tabOS[6,7] <- summary(m)$coef[3]
tabOS[6,8] <- summary(m)$conf.int[3]
tabOS[6,9] <- summary(m)$conf.int[4]
tabOS[6,10] <- summary(m)$coef[5]

# RFS

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,atrisk.title='Resection',legend.x="topright",legend.cex=1,legend.title='Resection',background=F,xlim=c(0,96),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))
title('RFS matching')

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

tabRFS[6,6] <- summary(m)$coef[2]
tabRFS[6,7] <- summary(m)$coef[3]
tabRFS[6,8] <- summary(m)$conf.int[3]
tabRFS[6,9] <- summary(m)$conf.int[4]
tabRFS[6,10] <- summary(m)$coef[5]

# semi-anatomic/wedge
data_sw<-data[data$PROCEDURE_cat=="semi-anatomic" | data$PROCEDURE_cat=="Wedge",]
data_sw$PROCEDURE_cat<-droplevels(data_sw$PROCEDURE_cat)
match <- matchit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                   INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                   MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                 data=data_sw,
                 method = "nearest", caliper=0.3)
data_complete_match <- match.data(match)

# Balance plots:

love.plot(match,threshold = .1,drop.distance = T,sample.names=c("Unmatched", "Matched"),
          stars="raw",stats="m",continuous="std",binary="std")

bal.tab(match, stats = "m", thresholds = .1, binary="std")

## MODELS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="bottomright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

##################################### BOOSTING ##############################################

# anatomic/semi-anatomic
data_tmp<-data_as
data_tmp$PROCEDURE_cat<-ifelse(data_tmp$PROCEDURE_cat=='semi-anatomic',1,0)
match <- matchit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                   INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                   MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                 data=data_tmp, distance="gbm", caliper = 0.3) # generalized boosted model
data_complete_match <- match.data(match)

# Balance plots:

love.plot(match,threshold = .1,drop.distance = T,sample.names=c("Unmatched", "Matched"),
          stars="raw",stats="m",continuous="std",binary="std")

bal.tab(match, stats = "m", thresholds = .1, binary="std")

## MODELS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="bottomright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

# anatomic/wedge
data_tmp<-data_aw
data_tmp$PROCEDURE_cat<-ifelse(data_tmp$PROCEDURE_cat=='Wedge',1,0)
match <- matchit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                   INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                   MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                 data=data_tmp,
                 method = "nearest", caliper=0.3)
data_complete_match <- match.data(match)

# Balance plots:

love.plot(match,threshold = .1,drop.distance = T,sample.names=c("Unmatched", "Matched"),
          stars="raw",stats="m",continuous="std",binary="std")

bal.tab(match, stats = "m", thresholds = .1, binary="std")

## MODELS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="bottomright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

# semi-anatomic/wedge
data_tmp<-data_sw
data_tmp$PROCEDURE_cat<-ifelse(data_tmp$PROCEDURE_cat=='Wedge',1,0)
match <- matchit(PROCEDURE_cat ~ eta_surg+bili_tot+Albuminemia+Creatininemia+
                   INR+Piastrine+Cirrosi+child_pugh_grade+Steatoepatite+inv_macro+
                   MVI+EdmondsonGrading+HBV+HCV+Potus.+satellitosi+R+sesso, 
                 data=data_sw,
                 method = "nearest", caliper=0.3)
data_complete_match <- match.data(match)

# Balance plots:

love.plot(match,threshold = .1,drop.distance = T,sample.names=c("Unmatched", "Matched"),
          stars="raw",stats="m",continuous="std",binary="std")

bal.tab(match, stats = "m", thresholds = .1, binary="std")

## MODELS

sFit <- prodlim(Hist(timeOSmesi,eventomorte)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Overall survival (%)",atrisk.cex=0.8,legend.x="bottomright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeOSmesi,eventomorte)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

sFit <- prodlim(Hist(timeRFSmesi,RFS)~PROCEDURE_cat,data=data_complete_match)
par(mar=c(3,1,1,1),mfrow=c(1,1))
plot(sFit,xlab="Time (Months)",ylab="Relapse Free survival (%)",atrisk.cex=0.8,legend.x="topright",legend.cex=1,background=F,logrank=T,xlim=c(0,60),axis1.at=seq(0,120,12))
summary(sFit,times = c(12,24,36,48,60,120))

m <- coxph(Surv(timeRFSmesi,RFS)~ PROCEDURE_cat+strata(subclass),data=data_complete_match)
summary(m)

#####################################################################################################

############################ Create summary tabs ########################################################

# create tables HR

write.csv(tabOS,'tabOS.csv',row.names = FALSE)
write.csv(tabRFS,'tabRFS.csv',row.names = FALSE)
