## THESIS ##
## Results simulations: plotting graphics

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

# SIMULATIONS: CASE A, B, C, D

# A
# read the table with the summary of the overall results
tabA<-read.csv('summary_A.csv')
row.names(tabA)<-c('adjusted','mult_weights','boost_weights','mult_gps','mult_classic')

# read the table with the overall results
tbA<-read.csv('sum_A.csv')

par(mfrow=c(2,2))
# BOXPLOT SEMI
HR<-c('s_HR_adj','s_HR_ipw_mn','s_HR_mn_gps','s_HR_mn_cl','s_HR_ipw_boost')
boxplot(tbA[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case A')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('s_HR_adj','s_HR_ipw_mn','s_HR_gps_mn','s_HR_cl_mn','s_HR_ipw_boost'), xpd = TRUE)
abline(h=0.8, lty=2, lwd=3)
mtext("Anatomic vs Semi-anatomic", side = 3, line = -1.5, outer = TRUE, font=2)

# BOXPLOT WEDGE
HR<-c('w_HR_adj','w_HR_ipw_mn','w_HR_mn_gps','w_HR_mn_cl','w_HR_ipw_boost')
boxplot(tbA[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case A')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('w_HR_adj','w_HR_ipw_mn','w_HR_gps_mn','w_HR_cl_mn','w_HR_ipw_boost'), xpd = TRUE)
abline(h=0.5, lty=2, lwd=3)
mtext("Anatomic vs Wedge", side = 3, line = -1.5, outer = TRUE, font=2)

################################################################################################
# B

# read the table with the summary of the overall results
tabB<-read.csv('summary_B.csv')
row.names(tabB)<-c('adjusted','mult_weights','boost_weights','mult_gps','mult_classic')

# read the table with the overall results
tbB<-read.csv('sum_B.csv')

# BOXPLOT sEMI
HR<-c('s_HR_adj','s_HR_ipw_mn','s_HR_mn_gps','s_HR_mn_cl','s_HR_ipw_boost')
boxplot(tbB[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case B')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('s_HR_adj','s_HR_ipw_mn','s_HR_gps_mn','s_HR_cl_mn','s_HR_ipw_boost'), xpd = TRUE)
abline(h=1.5, lty=2, lwd=3)

# BOXPLOT WEDGE
HR<-c('w_HR_adj','w_HR_ipw_mn','w_HR_mn_gps','w_HR_mn_cl','w_HR_ipw_boost')
boxplot(tbB[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case B')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('w_HR_adj','w_HR_ipw_mn','w_HR_gps_mn','w_HR_cl_mn','w_HR_ipw_boost'), xpd = TRUE)
abline(h=0.5, lty=2, lwd=3)

################################################################################################
# C
# read the table with the summary of the overall results
tabC<-read.csv('summary_C.csv')
row.names(tabC)<-c('adjusted','mult_weights','boost_weights','mult_gps','mult_classic')

# read the table with the overall results
tbC<-read.csv('sum_C.csv')

# BOXPLOT
HR<-c('s_HR_adj','s_HR_ipw_mn','s_HR_mn_gps','s_HR_mn_cl','s_HR_ipw_boost')
boxplot(tbC[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case C')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('s_HR_adj','s_HR_ipw_mn','s_HR_gps_mn','s_HR_cl_mn','s_HR_ipw_boost'), xpd = TRUE)
abline(h=1.5, lty=2, lwd=3)

# BOXPLOT WEDGE
HR<-c('w_HR_adj','w_HR_ipw_mn','w_HR_mn_gps','w_HR_mn_cl','w_HR_ipw_boost')
boxplot(tbC[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case C')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('w_HR_adj','w_HR_ipw_mn','w_HR_gps_mn','w_HR_cl_mn','w_HR_ipw_boost'), xpd = TRUE)
abline(h=2, lty=2, lwd=3)

################################################################################################
# D
# read the table with the summary of the overall results
tabD<-read.csv('summary_D.csv')
row.names(tabD)<-c('adjusted','mult_weights','boost_weights','mult_gps','mult_classic')

# read the table with the overall results
tbD<-read.csv('sum_D.csv')

# BOXPLOT
HR<-c('s_HR_adj','s_HR_ipw_mn','s_HR_mn_gps','s_HR_mn_cl','s_HR_ipw_boost')
boxplot(tbD[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case D')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('s_HR_adj','s_HR_ipw_mn','s_HR_gps_mn','s_HR_cl_mn','s_HR_ipw_boost'), xpd = TRUE)
abline(h=0.2, lty=2, lwd=3)

# BOXPLOT WEDGE
HR<-c('w_HR_adj','w_HR_ipw_mn','w_HR_mn_gps','w_HR_mn_cl','w_HR_ipw_boost')
boxplot(tbD[,HR],xaxt="n",xlab="", col = c('red','blue','purple','yellow','green'), main='Case D')
text(x =  seq_along(HR), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('w_HR_adj','w_HR_ipw_mn','w_HR_gps_mn','w_HR_cl_mn','w_HR_ipw_boost'), xpd = TRUE)
abline(h=0.7, lty=2, lwd=3)


############################################################################
# OVERALL SCATTERPLOT FOR THE ESTIMATIONS OF HRw and HRs

par(mfrow=c(1,1))
plot(c(1,2,3,4), c(tabA$HR_semi[1],tabB$HR_semi[1],tabC$HR_semi[1],tabD$HR_semi[1]), xaxt="n",  xlab="", ylab="", col=c('red'), type='b', bg='red', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_semi[2],tabB$HR_semi[2],tabC$HR_semi[2],tabD$HR_semi[2]), xaxt="n", xlab="", ylab="", col=c('blue'), type='b', bg='blue', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_semi[3],tabB$HR_semi[3],tabC$HR_semi[3],tabD$HR_semi[3]), xaxt="n", xlab="", ylab="", col=c('green'), type='b', bg='green', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_semi[4],tabB$HR_semi[4],tabC$HR_semi[4],tabD$HR_semi[4]), xaxt="n", xlab="", ylab="", col=c('purple'), bg='purple', type='b', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_semi[5],tabB$HR_semi[5],tabC$HR_semi[5],tabD$HR_semi[5]), xaxt="n", xlab="", ylab="", col=c('yellow'), bg='yellow', type='b', pch=16, lty=2)
abline(h=0.8, lty=2)
abline(h=1.5, lty=2)
abline(h=0.2, lty=2)
text(x =  seq_along(c('A','B','C','D')), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('A','B','C','D'), xpd = TRUE)
title('Mean of the estimated HRs')
legend(3.38, 1.6, legend=c("adj", "ipw_mn", "ipw_boost", "gps_mn", "cl_mn"),
       col=c("red", "blue", "green", "purple", "yellow"), pch=16)

par(mfrow=c(1,1))
plot(c(1,2,3,4), c(tabA$HR_wedged[1],tabB$HR_wedged[1],tabC$HR_wedged[1],tabD$HR_wedged[1]), xaxt="n", xlab="", ylab="", col=c('red'), type='b', bg='red', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_wedged[2],tabB$HR_wedged[2],tabC$HR_wedged[2],tabD$HR_wedged[2]), xaxt="n", xlab="", ylab="", col=c('blue'), type='b', bg='blue', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_wedged[3],tabB$HR_wedged[3],tabC$HR_wedged[3],tabD$HR_wedged[3]), xaxt="n", xlab="", ylab="", col=c('green'), type='b', bg='green', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_wedged[4],tabB$HR_wedged[4],tabC$HR_wedged[4],tabD$HR_wedged[4]), xaxt="n", xlab="", ylab="", col=c('purple'), bg='purple', type='b', pch=16, lty=2)
points(c(1,2,3,4), c(tabA$HR_wedged[5],tabB$HR_wedged[5],tabC$HR_wedged[5],tabD$HR_wedged[5]), xaxt="n", xlab="", ylab="", col=c('yellow'), bg='yellow', type='b', pch=16, lty=2)
abline(h=0.5, lty=2)
abline(h=2, lty=2)
abline(h=0.7, lty=2)
text(x =  seq_along(c('A','B','C','D')), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('A','B','C','D'), xpd = TRUE)
title('Mean of the estimated HRw')
legend(1, 1.9, legend=c("adj", "ipw_mn", "ipw_boost", "gps_mn", "cl_mn"),
       col=c("red", "blue", "green", "purple", "yellow"), pch=16)

##################################################################################################################
# LOLLIPLOTS FOR THE COVERAGE
# Unify all the tables
tab=rbind(tabA,tabB,tabC,tabD)

x_vec<-c(1,1.1,1.2,1.3,1.4,2,2.1,2.2,2.3,2.4,3,3.1,3.2,3.3,3.4,4,4.1,4.2,4.3,4.4)
col_cat<-c("red","blue", "green", "purple", "yellow","red","blue", "green", "purple", "yellow","red","blue", "green", "purple", "yellow","red","blue", "green", "purple", "yellow")
ggplot(tab, aes(x=x_vec, y=s_covarage)) +
        geom_segment( aes(x=x_vec, xend=x_vec, y=75, yend=s_covarage), color=col_cat) +
        geom_point( color=col_cat, size=4) +
        geom_hline(yintercept = c(95,100,75), linetype=c('dashed','solid','solid'), color = c('red','black','grey')) +
        geom_vline(xintercept = c(1.7,2.7,3.7), linetype='dashed', color = 'grey') +
        scale_x_continuous(breaks = c(1.2,2.2,3.2,4.2), labels = c('A','B','C','D')) +
        theme_light() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = 'top'
        ) +
        ylim(c(75,100)) +
        xlab("") +
        ylab("% of C.I. that contains the real value of HRs") +
        ggtitle('Coverage for semi-anatomic vs anatomic')

# wedge

ggplot(tab, aes(x=x_vec, y=w_covarage)) +
        geom_segment( aes(x=x_vec, xend=x_vec, y=75, yend=w_covarage), color=col_cat) +
        geom_point( color=col_cat, size=4) +
        geom_hline(yintercept = c(95,100,75), linetype=c('dashed','solid','solid'), color = c('red','black','grey')) +
        geom_vline(xintercept = c(1.7,2.7,3.7), linetype='dashed', color = 'grey') +
        scale_x_continuous(breaks = c(1.2,2.2,3.2,4.2), labels = c('A','B','C','D')) +
        theme_light() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = 'top'
        ) +
        ylim(c(75,100)) +
        xlab("") +
        ylab("% of C.I. that contains the real value of HRw") +
        ggtitle('Coverage for wedge vs anatomic')

##################################################################################################
############################# EMP STD ERROR VS MEAN SE LOG(HR) ######################################
# SCATTERPLOT

par(mfrow=c(1,1))
plot(c(1,2,3,4), c(1/(tabA$dev_semi/tabA$se_semi)[1],1/(tabB$dev_semi/tabB$se_semi)[1],1/(tabC$dev_semi/tabC$se_semi)[1],1/(tabD$dev_semi/tabD$se_semi)[1]), xaxt="n", xlab="", ylab="mean se of log(HR)/empirical se", ylim=c(0.65,1.1), col=c('red'), type='b', bg='red', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_semi/tabA$se_semi)[2],1/(tabB$dev_semi/tabB$se_semi)[2],1/(tabC$dev_semi/tabC$se_semi)[2],1/(tabD$dev_semi/tabD$se_semi)[2]), xaxt="n", xlab="", ylab="", col=c('blue'), type='b', bg='blue', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_semi/tabA$se_semi)[3],1/(tabB$dev_semi/tabB$se_semi)[3],1/(tabC$dev_semi/tabC$se_semi)[3],1/(tabD$dev_semi/tabD$se_semi)[3]), xaxt="n", xlab="", ylab="", col=c('green'), type='b', bg='green', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_semi/tabA$se_semi)[4],1/(tabB$dev_semi/tabB$se_semi)[4],1/(tabC$dev_semi/tabC$se_semi)[4],1/(tabD$dev_semi/tabD$se_semi)[4]), xaxt="n", xlab="", ylab="", col=c('purple'), bg='purple', type='b', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_semi/tabA$se_semi)[5],1/(tabB$dev_semi/tabB$se_semi)[5],1/(tabC$dev_semi/tabC$se_semi)[5],1/(tabD$dev_semi/tabD$se_semi)[5]), xaxt="n", xlab="", ylab="", col=c('yellow'), bg='yellow', type='b', pch=16, lty=2)
abline(h=1)
text(x =  seq_along(c('A','B','C','D')), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('A','B','C','D'), xpd = TRUE)
title('Mean std error of log(HRs) vs empirical std error', cex.main=0.9)
legend(3.2, 0.89, legend=c("adj", "ipw_mn", "ipw_boost", "gps_mn", "cl_mn"),
       col=c("red", "blue", "green", "purple", "yellow"), pch=16)

# wedge

par(mfrow=c(1,1))
plot(c(1,2,3,4), c(1/(tabA$dev_wedged/tabA$se_wedged)[1],1/(tabB$dev_wedged/tabB$se_wedged)[1],1/(tabC$dev_wedged/tabC$se_wedged)[1],1/(tabD$dev_wedged/tabD$se_wedged)[1]), xaxt="n", xlab="", ylab="mean se of log(HR)/empirical se", ylim=c(0.6,1.1), col=c('red'), type='b', bg='red', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_wedged/tabA$se_wedged)[2],1/(tabB$dev_wedged/tabB$se_wedged)[2],1/(tabC$dev_wedged/tabC$se_wedged)[2],1/(tabD$dev_wedged/tabD$se_wedged)[2]), xaxt="n", xlab="", ylab="", col=c('blue'), type='b', bg='blue', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_wedged/tabA$se_wedged)[3],1/(tabB$dev_wedged/tabB$se_wedged)[3],1/(tabC$dev_wedged/tabC$se_wedged)[3],1/(tabD$dev_wedged/tabD$se_wedged)[3]), xaxt="n", xlab="", ylab="", col=c('green'), type='b', bg='green', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_wedged/tabA$se_wedged)[4],1/(tabB$dev_wedged/tabB$se_wedged)[4],1/(tabC$dev_wedged/tabC$se_wedged)[4],1/(tabD$dev_wedged/tabD$se_wedged)[4]), xaxt="n", xlab="", ylab="", col=c('purple'), bg='purple', type='b', pch=16, lty=2)
points(c(1,2,3,4), c(1/(tabA$dev_wedged/tabA$se_wedged)[5],1/(tabB$dev_wedged/tabB$se_wedged)[5],1/(tabC$dev_wedged/tabC$se_wedged)[5],1/(tabD$dev_wedged/tabD$se_wedged)[5]), xaxt="n", xlab="", ylab="", col=c('yellow'), bg='yellow', type='b', pch=16, lty=2)
abline(h=1)
text(x =  seq_along(c('A','B','C','D')), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,
     labels = c('A','B','C','D'), xpd = TRUE)
title('Mean std error of log(HRs) vs empirical std error', cex.main=0.9)
legend(1, 0.9, legend=c("adj", "ipw_mn", "ipw_boost", "gps_mn", "cl_mn"),
       col=c("red", "blue", "green", "purple", "yellow"), pch=16)

####################################################################################################
########### SUMMARY TABLES ##########################################################################

##### A

ts<-tabA[,c(1,2,3,7,8)]
ts$ratio_semi<-ts$se_semi/ts$dev_semi
ts$s_coverage<-tabA$s_covarage

tw<-tabA[,c(4,5,6,10,11)]
tw$ratio_wedged<-tw$se_wedged/tw$dev_wedged
tw$w_coverage<-tabA$w_covarage

t<-ts
t[,c(8:14)]<-tw
write.csv(t,file='tableA.csv')

##### B

ts<-tabB[,c(1,2,3,7,8)]
ts$ratio_semi<-ts$se_semi/ts$dev_semi
ts$s_coverage<-tabB$s_covarage

tw<-tabB[,c(4,5,6,10,11)]
tw$ratio_wedged<-tw$se_wedged/tw$dev_wedged
tw$w_coverage<-tabB$w_covarage

t<-ts
t[,c(8:14)]<-tw
write.csv(t,file='tableB.csv')

##### C

ts<-tabC[,c(1,2,3,7,8)]
ts$ratio_semi<-ts$se_semi/ts$dev_semi
ts$s_coverage<-tabC$s_covarage

tw<-tabC[,c(4,5,6,10,11)]
tw$ratio_wedged<-tw$se_wedged/tw$dev_wedged
tw$w_coverage<-tabC$w_covarage

t<-ts
t[,c(8:14)]<-tw
write.csv(t,file='tableC.csv')

##### D

ts<-tabD[,c(1,2,3,7,8)]
ts$ratio_semi<-ts$se_semi/ts$dev_semi
ts$s_coverage<-tabD$s_covarage

tw<-tabD[,c(4,5,6,10,11)]
tw$ratio_wedged<-tw$se_wedged/tw$dev_wedged
tw$w_coverage<-tabD$w_covarage

t<-ts
t[,c(8:14)]<-tw
write.csv(t,file='tableD.csv')
