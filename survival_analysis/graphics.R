## THESIS ##
## Graphics SUrvival Analysis ##

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

# set the path
setwd("name_of_the_path")

##################################### FOREST PLOT ##############################################

# create table with overall results
tb<-as.data.frame(matrix(nrow=5,ncol=6))
row.names(tb)<-c('adjusted','mult_weights','boost_weights','mult_gps','boost_gps','mult_classic', 'boost_classic')
colnames(tb)<-c('HR_semi','lower_semi','upper_semi','HR_wedge','lower_wedge','upper_wedge')

HR_semi<-c()
lower_semi<-c()
upper_semi<-c()
HR_wedged<-c()
lower_wedged<-c()
upper_wedged<-c()

# compute the mean of the values obtained
tab_pre<-read.csv('sum.csv')
for(j in 1:ncol(tab_pre)){
  if(grepl('s_HR',colnames(tab_pre)[j])){
    HR_semi<-c(HR_semi,mean(tab_pre[,j],na.rm=TRUE))
  }
  if(grepl('w_HR',colnames(tab_pre)[j])){
    HR_wedged<-c(HR_wedged,mean(tab_pre[,j],na.rm=TRUE))
  }
  if(grepl('s_LOWER',colnames(tab_pre)[j])){
    lower_semi<-c(lower_semi,mean(tab_pre[,j],na.rm=TRUE))
  }
  if(grepl('w_LOWER',colnames(tab_pre)[j])){
    lower_wedged<-c(lower_wedged,mean(tab_pre[,j],na.rm=TRUE))
  }
  if(grepl('s_UPPER',colnames(tab_pre)[j])){
    upper_semi<-c(upper_semi,mean(tab_pre[,j],na.rm=TRUE))
  }
  if(grepl('w_UPPER',colnames(tab_pre)[j])){
    upper_wedged<-c(upper_wedged,mean(tab_pre[,j],na.rm=TRUE))
  }
}

tb[,1]<-HR_semi
tb[,2]<-lower_semi
tb[,3]<-upper_semi
tb[,4]<-HR_wedged
tb[,5]<-lower_wedged
tb[,6]<-upper_wedged

### FOREST PLOTS ###

# OS
fp <- ggplot(data=tb[,c(1,2,3)], aes(x=rownames(tb), y=HR_semi, ymin=lower_semi, ymax=upper_semi)) +
  geom_pointrange() + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("HR (95% CI)") +
  theme_bw()  # use a white background
print(fp)

# RFS
fp <- ggplot(data=tb[,c(4,5,6)], aes(x=rownames(tb), y=HR_wedged, ymin=lower_wedged, ymax=upper_wedged)) +
  geom_pointrange() + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("HR (95% CI)") +
  theme_bw()  # use a white background
print(fp)

# dev std beta (empirical std error) vs mean std error log(HR)
dev_semi<-c()
dev_wedged<-c()
se_semi<-c()
se_wedged<-c()

for(j in 1:ncol(tab_pre)){
  if(grepl('s_HR',colnames(tab_pre)[j])){
    dev_semi<-c(dev_semi,sd(log(tab_pre[,j]),na.rm=TRUE))
  }
  if(grepl('w_HR',colnames(tab_pre)[j])){
    dev_wedged<-c(dev_wedged,sd(log(tab_pre[,j]),na.rm=TRUE))
  }
  if(grepl('s_SE',colnames(tab_pre)[j])){
    se_semi<-c(se_semi,mean(tab_pre[,j],na.rm=TRUE))
  }
  if(grepl('w_SE',colnames(tab_pre)[j])){
    se_wedged<-c(se_wedged,mean(tab_pre[,j],na.rm=TRUE))
  }
}

tb$dev_semi<-dev_semi
tb$se_semi<-se_semi
tb$diff_semi<-abs(tb$dev_semi-tb$se_semi)
tb$dev_wedged<-dev_wedged
tb$se_wedged<-se_wedged
tb$diff_wedged<-abs(tb$dev_wedged-tb$se_wedged)

tab_pre<-tab_pre[!is.na(tab_pre$s_UPPER_adj),]
HRs=0.2
HRw=0.7
s_covarage<-c()
s_covarage<-c(s_covarage,nrow(tab_pre[tab_pre['s_LOWER_adj']<= HRs & HRs <= tab_pre['s_UPPER_adj'],])/nrow(tab_pre)*100)
s_covarage<-c(s_covarage,nrow(tab_pre[tab_pre['s_LOWER_ipw_mn']<= HRs & HRs <= tab_pre['s_UPPER_ipw_mn'],])/nrow(tab_pre)*100)
s_covarage<-c(s_covarage,nrow(tab_pre[tab_pre['s_LOWER_ipw_boost']<= HRs & HRs <= tab_pre['s_UPPER_ipw_boost'],])/nrow(tab_pre)*100)
s_covarage<-c(s_covarage,nrow(tab_pre[tab_pre['s_LOWER_mn_gps']<= HRs & HRs <= tab_pre['s_UPPER_mn_gps'],])/nrow(tab_pre)*100)
s_covarage<-c(s_covarage,nrow(tab_pre[tab_pre['s_LOWER_mn_cl']<= HRs & HRs <= tab_pre['s_UPPER_mn_cl'],])/nrow(tab_pre)*100)
w_covarage<-c()
w_covarage<-c(w_covarage,nrow(tab_pre[tab_pre['w_LOWER_adj']<= HRw & HRw <= tab_pre['w_UPPER_adj'],])/nrow(tab_pre)*100)
w_covarage<-c(w_covarage,nrow(tab_pre[tab_pre['w_LOWER_ipw_mn']<= HRw & HRw <= tab_pre['w_UPPER_ipw_mn'],])/nrow(tab_pre)*100)
w_covarage<-c(w_covarage,nrow(tab_pre[tab_pre['w_LOWER_ipw_boost']<= HRw & HRw <= tab_pre['w_UPPER_ipw_boost'],])/nrow(tab_pre)*100)
w_covarage<-c(w_covarage,nrow(tab_pre[tab_pre['w_LOWER_mn_gps']<= HRw & HRw <= tab_pre['w_UPPER_mn_gps'],])/nrow(tab_pre)*100)
w_covarage<-c(w_covarage,nrow(tab_pre[tab_pre['w_LOWER_mn_cl']<= HRw & HRw <= tab_pre['w_UPPER_mn_cl'],])/nrow(tab_pre)*100)

tb$s_covarage<-s_covarage
tb$w_covarage<-w_covarage

write.csv(tb, file="summary.csv", row.names = FALSE)
