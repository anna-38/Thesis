## THESIS ##
## descriptive analysis ##
## cleaning data ##

rm(list=ls())

# import libraries
library(xlsx)
library(ggplot2)
library(MatchIt)
library(cobalt)
library(prodlim)
library(survival)
library(mice) # for multiple imputation

# set the path
setwd("name_of_the_path")

# read the dataset
data <-read.xlsx("large HCC.xlsx", sheetName = "large_HCC")
str(data)
View(data)

print(paste('Numero osservazioni =', length(data$ID_Paziente), sep=' '))
print(paste('Numero variabili =', length(data), sep=' '))

### MODIFY DATASET ###
# transform the categorical variables into factor
var<-c('ID_CENTRO','ID_Paziente','sesso','ASA','CCI','Cirrosi','Steatoepatite','child_pugh_grade','MELD_score','HBV','HCV','eradicato_con_DDA','Potus.','inv_macro_port','inv_macro_sovraep','BCLC.','neoadj','ARvsPSR','PROCEDURE','Laparoscopia','conversione','eco_intraop','trombosi_portale','Pringle','bloodlos','trasfusione_intraopUEC','sopravvissuto_chirurgia','complicanze','clavien_dindo','ascite_postop','X90_d_mortality','reop_complicanze','istotipo','EdmondsonGrading','MVI','satellitosi','capsula','R','causadecesso','eventorecidiva','single_multipleREC','bilobarREC','localizzazioneREC','localREC','nnoduliREC','sizeREC','AFPREC','Rec_treat','X..DELTA','size.radio...35.mm','MAX...49')
# RFS and eventomorte are NOT included, otherwise they'll give problem with Cox model
data[,var]<-lapply(data[,var], function(x) factor(x))
# eventomorte is character -> convert to numeric
data$eventomorte<-as.numeric(data$eventomorte)
data$Alfa.Feto.Protein<-as.numeric(data$Alfa.Feto.Protein)
data$bloodlos<-as.numeric(data$bloodlos)
summary(data)

à####### to create descriptive statistics table for excel ##############
library(tableone)
table  <- CreateTableOne(vars =
                           names(subset(data[4:length(data)],select=-c(data_surgery,data_dimissione,X90_d_mortality,dataMorte_ultimoFU,datarecidiva_ultimoFU))),
                         data=data,includeNA=F)
summary(table)
print(table,quote=T,missing=T)
#################################################

################## DESCRIPTIVE ANALYSIS OF CORE VARIABLES: ##############################
########################## OUTPUTS AND PROCEDURE ########################################

# OUTCOMES: RFS timeRFSmesi eventomorte timeOSmesi

# RFS
summary(as.factor(data$RFS))

df<-as.data.frame(summary(as.factor(data[!is.na(data$RFS),]$RFS))/length(data[!is.na(data$RFS),]$RFS)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(75,25), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to RFS') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# OS
summary(as.factor(data$eventomorte))

df<-as.data.frame(summary(as.factor(data[!is.na(data$eventomorte),]$eventomorte))/length(data[!is.na(data$eventomorte),]$eventomorte)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(70,20), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to eventomorte') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# timeRFSmesi
summary(data$timeRFSmesi)
# Regisetered times (in years)
boxplot(data$timeRFSmesi/12, data$timeRFSmesi[data$RFS==0]/12, data$timeRFSmesi[data$RFS==1]/12, ylab="time (in years)", names=c('all','censored','events'), col=c('orange','violet','lightblue'), main="Time RFS in years")

# timeOSmesi
summary(data$timeOSmesi)
# Registered times (in years)
boxplot(data$timeOSmesi/12, data$timeOSmesi[data$eventomorte==0]/12, data$timeOSmesi[data$eventomorte==1]/12, ylab="time (in mesi)", names=c('all','censored','events'), col=c('orange','violet','lightblue'), main="Time OS in years")
# timeRFSmesi==0 check:
nrow(data[data$timeRFSmesi==0 & !is.na(data$timeRFSmesi),]) # only 1
data[data$timeRFSmesi==0 & !is.na(data$timeRFSmesi),] # who?
# modify the 0s
data$timeRFSmesi[data$timeRFSmesi==0]=1/30
data$timeOSmesi[data$timeOSmesi==0]=1/30

# PROCEDURE
summary(data$PROCEDURE)

df<-as.data.frame(summary(data[!is.na(data$PROCEDURE),]$PROCEDURE)/length(data[!is.na(data$PROCEDURE),]$PROCEDURE)*100)
# Pie chart
PROCEDURE<-rownames(df)
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=PROCEDURE)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  ggtitle('% patients divided in respect to procedure') +
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

proc<-data.frame(matrix(nrow=6,ncol=1), row.names=rownames(df))
colnames(proc)='#patients'
proc[,1]<-df[,1]*length(data[-c(as.numeric(rownames(data[is.na(data$PROCEDURE),]))),]$PROCEDURE)/100
# 100% stacked bar chart with groups of age
ggplot() + aes(proc, fill=rownames(proc), y = colnames(proc), x = proc[,1]/sum(proc[,1])*100) +
  geom_bar(stat="identity", position=position_stack(reverse = TRUE)) +
  guides(fill=guide_legend(title="")) +
  scale_fill_brewer(palette="Set2") +
  ggtitle("Patients grouped by procedure") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") +
  xlab("% patients")
# n patients/percentage patients
proc
proc/sum(proc)*100

# create with three categories
data$PROCEDURE_cat<-as.factor(ifelse(data$PROCEDURE==0,"Wedge",ifelse(data$PROCEDURE==4,"anatomic","semi-anatomic")))

df<-as.data.frame(summary(data[!is.na(data$PROCEDURE_cat),]$PROCEDURE_cat)/length(data[!is.na(data$PROCEDURE_cat),]$PROCEDURE_cat)*100)
proc<-data.frame(matrix(nrow=3,ncol=1), row.names=rownames(df))
colnames(proc)='#patients'
proc[,1]<-df[,1]*length(data[-c(as.numeric(rownames(data[is.na(data$PROCEDURE_cat),]))),]$PROCEDURE_cat)/100
# 100% stacked bar chart with groups of age
ggplot() + aes(proc, fill=rownames(proc), y = colnames(proc), x = proc[,1]/sum(proc[,1])*100) +
  geom_bar(stat="identity", position=position_stack(reverse = TRUE)) +
  guides(fill=guide_legend(title="")) +
  scale_fill_brewer(palette="Set2") +
  ggtitle("Patients grouped by procedure") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") +
  xlab("% patients")
# n patients/percentage patients
proc
proc/sum(proc)*100

############################# MODIFY DATASET #################################################

# REMOVE:
# - PROCEDURE==NA: only one, obviously it's a problem
# - PROCEDURE==5: They are "other treatments", we are not interested in them
# - RFS==NA/eventomorte==NA: we need the outcome
# - timeRFSmesi/timeOSmesi==NA: we need the outcome

data_complete<-data # 613 instances

# -PROCEDURE==NA
nrow(data[is.na(data$PROCEDURE),]) # how many instances will be eliminated? 1
data<-data[-c(as.numeric(rownames(data[is.na(data$PROCEDURE),]))),] #612 rows

# -PROCEDURE==5
nrow(data[data$PROCEDURE==5,])
data<-data[data$PROCEDURE!=5,]
data$PROCEDURE<-droplevels(data$PROCEDURE)

# - RFS/eventomorte==NA
nrow(data[is.na(data$RFS) | is.na(data$eventomorte),]) # how many instances will be eliminated? 13
data<-data[!is.na(data$RFS) & !is.na(data$eventomorte),] #555 rows

# - timeRFSmesi/timeOSmesi==NA
nrow(data[is.na(data$timeOSmesi) | is.na(data$timeRFSmesi),]) # how many instances will be eliminated? 17
data<-data[!is.na(data$timeOSmesi) & !is.na(data$timeRFSmesi),] #538 rows

nrow(data) # 538 rows (vs starting 613 rows)

####################### DESCRIPTIVE ANALYSIS ################################################

# AGE
summary(data$eta_surg)

# <18 SHOULDN'T be part of the dataset! --> remove them
data<-data[data$eta_surg>=18 | is.na(data$eta_surg),]
nrow(data) # final dataset: 536 rows (vs starting 613 rows)

# Divide in groups for age
gruppi<-factor(c('age<40','40<=age<50','50<=age<60','60<=age<70',"70<=age<80","80=<age"),levels=c('age<40','40<=age<50','50<=age<60','60<=age<70',"70<=age<80","80=<age"))
eta<-data.frame(matrix(nrow=6,ncol=1), row.names=gruppi)
colnames(eta)='number patients'
# exclude the row with NA
val<-c(sum(as.numeric(data[!is.na(data$eta_surg),]$eta_surg<40)),sum(as.numeric(40<=data[!is.na(data$eta_surg),]$eta_surg & data[!is.na(data$eta_surg),]$eta_surg<50)),sum(as.numeric(50<=data[!is.na(data$eta_surg),]$eta_surg & data[!is.na(data$eta_surg),]$eta_surg<60)),sum(as.numeric(60<=data[!is.na(data$eta_surg),]$eta_surg & data[!is.na(data$eta_surg),]$eta_surg<70)),sum(as.numeric(70<=data[!is.na(data$eta_surg),]$eta_surg & data[!is.na(data$eta_surg),]$eta_surg<80)),sum(as.numeric(80<=data[!is.na(data$eta_surg),]$eta_surg)))
eta[,1]<-val
# 100% stacked bar chart with groups of age
ggplot() + aes(eta, fill=gruppi, y = colnames(eta), x = val/sum(val)*100) +
  geom_bar(stat="identity", position=position_stack(reverse = TRUE)) +
  guides(fill=guide_legend(title="")) +
  scale_fill_brewer(palette="Set2") +
  ggtitle("Patients grouped by age") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") +
  xlab("% patients")
# n patients/percentage patients
eta
eta/sum(eta)*100

# Boxplot
boxplot(data$eta_surg, main="Age of patients")

# SEX
# Sex: absolute and relative frequencies
summary(data$sesso)
df<-as.data.frame(summary(data$sesso)/length(data$sesso)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(90,40), label = rownames(df)), color = "black", size=4.5) +
  ggtitle('% patients divided in respect to sex') +
  scale_fill_manual(values=c("pink","lightblue")) +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.4),
        # plot.title = element_text(face = "bold"),
        panel.background = element_blank()) +
  xlab("") +
  ylab("")

# ASA - many NA
summary(data$ASA)

# remove NA
df<-as.data.frame(summary(data[!is.na(data$ASA),]$ASA)/length(data[!is.na(data$ASA),]$ASA)*100)
# Pie chart
ASA<-rownames(df)
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=ASA)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  ggtitle('% patients divided in respect to ASA') +
  theme(legend.position="right") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_brewer(palette="Set2") +
  xlab("") +
  ylab("")

boxplot(as.numeric(data$ASA))

# CCI - many NA
summary(data$CCI)
# since there aren't instances with CCI=0 anymore:
data$CCI<-droplevels(data$CCI)
df<-as.data.frame(summary(data[!is.na(data$CCI),]$CCI)/length(data[!is.na(data$CCI),]$CCI)*100)
df

boxplot(as.numeric(data$CCI), main="CCI")

# MELD score - many NA
summary(data$MELD_score)
boxplot(as.numeric(data$MELD_score), main="MELD Score")

# cirrosi
summary(data$Cirrosi)

df<-as.data.frame(summary(data[!is.na(data$Cirrosi),]$Cirrosi)/length(data[!is.na(data$Cirrosi),]$Cirrosi)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(75,20), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to Cirrosi') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# Child purgh grade
summary(data$child_pugh_grade)

df<-as.data.frame(summary(data[!is.na(data$child_pugh_grade),]$child_pugh_grade)/length(data[!is.na(data$child_pugh_grade),]$child_pugh_grade)*100)
# Pie chart
Child<-rownames(df)
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=Child)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  ggtitle('% patients divided in respect to Child') +
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# unify 1 and 2
data[data$child_pugh_grade==2 & !is.na(data$child_pugh_grade),]$child_pugh_grade<-1
data$child_pugh_grade<-droplevels(data$child_pugh_grade)

# Steatoepatite
summary(data$Steatoepatite)

df<-as.data.frame(summary(data[!is.na(data$Steatoepatite),]$Steatoepatite)/length(data[!is.na(data$Steatoepatite),]$Steatoepatite)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(70,15), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to Steatoepatite') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# HCV
summary(data$HCV)

df<-as.data.frame(summary(data[!is.na(data$HCV),]$HCV)/length(data[!is.na(data$HCV),]$HCV)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(70,15), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to HCV') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# HBV
summary(data$HBV)

df<-as.data.frame(summary(data[!is.na(data$HBV),]$HBV)/length(data[!is.na(data$HBV),]$HBV)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(60,8), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to HBV') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# potus
summary(data$Potus.)

df<-as.data.frame(summary(data[!is.na(data$Potus.),]$Potus.)/length(data[!is.na(data$Potus.),]$Potus.)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(60,8), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to potus') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# bilirubina
summary(data$bili_tot)
boxplot(data$bili_tot, main="Bilirubina")

# albuminemia
summary(data$Albuminemia)
boxplot(data$Albuminemia, main="Albuminemia")

# creatinina
summary(data$Creatininemia)
boxplot(data$Creatininemia, main="Creatinina")

# piastrine
summary(data$Piastrine)
boxplot(data$Piastrine, main="Piastrine")

# INR
summary(data$INR)
boxplot(data$INR, main="INR")

# inv_macro_port
summary(data$inv_macro_port)

df<-as.data.frame(summary(data[!is.na(data$inv_macro_port),]$inv_macro_port)/length(data[!is.na(data$inv_macro_port),]$inv_macro_port)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(55,7), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to inv_macro_port') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# inv_macro_sovraep
summary(data$inv_macro_sovraep)

df<-as.data.frame(summary(data[!is.na(data$inv_macro_sovraep),]$inv_macro_sovraep)/length(data[!is.na(data$inv_macro_sovraep),]$inv_macro_sovraep)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(55,5), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to inv_macro_sovraep') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# EdmondsonGrading
summary(data$EdmondsonGrading)

# remove NA
df<-as.data.frame(summary(data[!is.na(data$EdmondsonGrading),]$EdmondsonGrading)/length(data[!is.na(data$EdmondsonGrading),]$EdmondsonGrading)*100)
# Pie chart
EdmondsonGrading=rownames(df)
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=EdmondsonGrading)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  ggtitle('% patients divided in respect to EdmondsonGrading') +
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# MVI - many NA
summary(data$MVI)
# Problem: we have white spaces, they need to be substituted by NA
data[data$MVI==" " & !is.na(data$MVI),]$MVI<-NA

# remove NA
df<-as.data.frame(summary(data[!is.na(data$MVI),]$MVI)/length(data[!is.na(data$MVI),]$MVI)*100)
# Pie chart
MVI<-rownames(df)[-1]
ggplot(as.data.frame(df[-1,]), aes(x="", y=df[-1,1], fill=MVI)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  ggtitle('% pazienti divisi rispetto a MVI') +
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# unify 1 and 2
data[data$MVI==2 & !is.na(data$MVI),]$MVI<-1
data$MVI<-droplevels(data$MVI)

# satellitosi
summary(data$satellitosi)

# remove NA
df<-as.data.frame(summary(data[!is.na(data$satellitosi),]$satellitosi)/length(data[!is.na(data$satellitosi),]$satellitosi)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(60,10), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to satellitosi') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# capsula
summary(data$capsula)

# remove NA
df<-as.data.frame(summary(data[!is.na(data$capsula),]$capsula)/length(data[!is.na(data$capsula),]$capsula)*100)
# Pie chart
ggplot(as.data.frame(df), aes(x="", y=df[,1], fill=rownames(df))) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = c(70,20), label = rownames(df)), color = "white", size=4.5) +
  ggtitle('% patients divided in respect to capsula') +
  theme(legend.position="none") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# R - many NA
summary(data$R)

# Problem: we have white spaces, they need to be substituted by NA
data[data$R==" " & !is.na(data$R),]$R<-NA

# remove NA
df<-as.data.frame(summary(data[!is.na(data$R),]$R)/length(data[!is.na(data$R),]$R)*100)
# Pie chart
R<-rownames(df)[-1]
ggplot(as.data.frame(df[-1,]), aes(x="", y=df[-1,1], fill=R)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  ggtitle('% pazienti divisi rispetto a R') +
  theme(legend.position="right") +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  xlab("") +
  ylab("")

# unify 1 and 2
data[data$R==2 & !is.na(data$R),]$R<-1
data$R<-droplevels(data$R)

############## UNIVARIATE ANALYSIS #######################
# To check which variables seem to have the strongest correlation with the output
vars<-c("eta_surg",
        "sesso",
        "ASA",
        "CCI",
        "MELD_score",
        "Cirrosi",
        "child_pugh_grade",
        "Steatoepatite",
        "HCV",
        "HBV",
        "Potus.",
        "bili_tot",
        "Albuminemia",
        "Creatininemia",
        "Piastrine",
        "INR",
        "inv_macro_port",
        "inv_macro_sovraep",
        "PROCEDURE",
        "EdmondsonGrading",
        "MVI",
        "satellitosi",
        "capsula",
        "R")

# check if there are constant values:
tab<-ifelse(n <- sapply(vars, function(x) length(levels(data[,x]))) == 1, "DROP", "NODROP")
print('Variables that are constant - they will be removed: ')
tab[tab=="DROP"]

results<-as.data.frame(matrix(nrow=length(vars),ncol=5,byrow=T))
k=1
for (var in vars) {
  formula    <- as.formula(paste("Surv(timeRFSmesi,RFS)~",var))
  fit.uni <- coxph(formula,data)
  beta <- coef(fit.uni)
  se   <- sqrt(diag(fit.uni$var))
  CI   <- round(exp(confint(fit.uni)), 3)
  features <- names(fit.uni$coefficients)
  if (length(beta)==1) {
    results[k,]=c(round(c(exp(beta), CI, p=1-pchisq((beta/se)^2, 1)),3),features)
    k=k+1
  } else {
    for (i in range(1,length(beta))) {
      results[k,]=c(round(c(exp(beta[i]), CI[i,], p=1-pchisq((beta[i]/se[i])^2, 1)),3),features[i])
      k=k+1
    } 
  }
}
names(results)<-c("HR","lower95%CI","upper95%CI","p","feature") # HR=hazard ratio
print(results)

write.csv(results, 'univariate_cox.csv', row.names = FALSE)

# Which var significant in the univariate analysis?
print("Variables which result to be significative in the univariate analysis: ")
print(results$feature[results$p<0.05])
sig_vars<-results$feature[results$p<0.05]

############################### MISSING VALUES: MULTIPLE IMPUTATION ###############################################
# Fill NA values with values computed with multiple imputation:
# NB even the OUTPUT variables must be put in the model.
# If there is a connection, even variables NOT in the final model!!
# variables selected for the model
vars_mi<-c("eta_surg",
           "sesso",
           "Cirrosi",
           "child_pugh_grade",
           "Steatoepatite",
           "HCV",
           "HBV",
           "Potus.",
           "bili_tot",
           "Albuminemia",
           "Creatininemia",
           "INR",
           "Piastrine",
           "inv_macro_port",
           "inv_macro_sovraep",
           "MVI",
           "EdmondsonGrading",
           "satellitosi",
           "R",
           "RFS",
           "eventomorte",
           "timeRFSmesi",
           "timeOSmesi",
           "PROCEDURE_cat")
vars_mi

datami<-data[,vars_mi]

# Pattern of missingness:
# there is some subjects a bit problematic, with many NA
# (NA is indicated by purple squares)
par(mfrow=c(1,1),mar=c(1,1,1,1))
md.pattern(datami)

my.m<-5 # 5 is the default number of multiple imputations

imp <- mice(datami,m=my.m,maxit=0,seed=123, print=FALSE)
imp <- mice(datami,m=my.m,maxit=21,seed=123, print=FALSE)

# averaging imputed datasets
# stack data
data.long<-complete(imp,"long")

#### Save imputed dataset: not to compute imputation everytime ####
# write.table(data.long,"large HCC long.txt")

# Read the saved dataset
data.long<-read.table("large HCC long.txt",header=T)


# take the mean for quantitative vars and the mode for categorical vars

# functions to find the mode (if many modes, pick one randomly)
my.mode <- function(x) {
  mode<-names(table(x)[table(x)==max(table(x))])
  mode<-ifelse(length(mode)==1,mode,sample(mode,1))
  if(length(mode)>1){
    print("WARNING: more than one mode, mode picked randomly")
  }
  return(mode)
}

# create new dataset
data.avg<-data.frame(id=1:(dim(data.long)[1]/5))

# numeric
data.avg$eta_surg <- as.numeric(with(data.long, by(eta_surg, as.factor(.id), FUN=mean)))
data.avg$bili_tot <- as.numeric(with(data.long, by(bili_tot, as.factor(.id), FUN=mean)))
data.avg$Albuminemia <- as.numeric(with(data.long, by(Albuminemia, as.factor(.id), FUN=mean)))
data.avg$Creatininemia <- as.numeric(with(data.long, by(Creatininemia, as.factor(.id), FUN=mean)))
data.avg$INR <- as.numeric(with(data.long, by(INR, as.factor(.id), FUN=mean)))
data.avg$Piastrine <- as.numeric(with(data.long, by(Piastrine, as.factor(.id), FUN=mean)))

# categorical
data.avg$Cirrosi <- as.factor(with(data.long, by(Cirrosi, as.factor(.id), FUN=my.mode)))
data.avg$child_pugh_grade <- as.factor(with(data.long, by(child_pugh_grade, as.factor(.id), FUN=my.mode)))
data.avg$Steatoepatite <- as.factor(with(data.long, by(Steatoepatite, as.factor(.id), FUN=my.mode)))
data.avg$inv_macro_port <- as.factor(with(data.long, by(inv_macro_port, as.factor(.id), FUN=my.mode)))
data.avg$inv_macro_sovraep <- as.factor(with(data.long, by(inv_macro_sovraep, as.factor(.id), FUN=my.mode)))
data.avg$MVI <- as.factor(with(data.long, by(MVI, as.factor(.id), FUN=my.mode)))
data.avg$EdmondsonGrading <- as.factor(with(data.long, by(EdmondsonGrading, as.factor(.id), FUN=my.mode)))
data.avg$HBV <- as.factor(with(data.long, by(HBV, as.factor(.id), FUN=my.mode)))
data.avg$HCV <- as.factor(with(data.long, by(HCV, as.factor(.id), FUN=my.mode)))
data.avg$Potus. <- as.factor(with(data.long, by(Potus., as.factor(.id), FUN=my.mode)))
data.avg$satellitosi <- as.factor(with(data.long, by(satellitosi, as.factor(.id), FUN=my.mode)))
data.avg$R<-as.factor(with(data.long, by(R, as.factor(.id), FUN=my.mode)))

# remained unmodified
data.avg$sesso<-data$sesso
data.avg$PROCEDURE<-data$PROCEDURE
data.avg$PROCEDURE_cat<-data$PROCEDURE_cat
data.avg$RFS<-data$RFS
data.avg$timeRFSmesi<-data$timeRFSmesi
data.avg$eventomorte<-data$eventomorte
data.avg$timeOSmesi<-data$timeOSmesi

# data.avg NEW DATASET

write.csv(data.avg, file="data.avg.csv", row.names = FALSE)

############################ UNIVARIATE ANALYSIS WITH COMPLETE DATASET ####################

vars_ua<-vars_mi<-c("eta_surg",
                    "sesso",
                    "Cirrosi",
                    "child_pugh_grade",
                    "Steatoepatite",
                    "HCV",
                    "HBV",
                    "Potus.",
                    "bili_tot",
                    "Albuminemia",
                    "Creatininemia",
                    "INR",
                    "Piastrine",
                    "inv_macro_port",
                    "inv_macro_sovraep",
                    "MVI",
                    "EdmondsonGrading",
                    "satellitosi",
                    "R",
                    "PROCEDURE")

# check if there are constant values:
tab<-ifelse(n <- sapply(vars_ua, function(x) length(levels(data.avg[,x]))) == 1, "DROP", "NODROP")
print('Variables that are constant - they will be removed: ')
tab[tab=="DROP"] # there arent't

results<-as.data.frame(matrix(nrow=length(vars_ua),ncol=5,byrow=T))
k=1
for (var in vars_ua) {
  formula    <- as.formula(paste("Surv(timeRFSmesi,RFS)~",var))
  fit.uni <- coxph(formula,data.avg)
  beta <- coef(fit.uni)
  se   <- sqrt(diag(fit.uni$var))
  CI   <- round(exp(confint(fit.uni)), 3)
  features <- names(fit.uni$coefficients)
  if (length(beta)==1) {
    results[k,]=c(round(c(exp(beta), CI, p=1-pchisq((beta/se)^2, 1)),3),features)
    k=k+1
  } else {
    for (i in range(1,length(beta))) {
      results[k,]=c(round(c(exp(beta[i]), CI[i,], p=1-pchisq((beta[i]/se[i])^2, 1)),3),features[i])
      k=k+1
    } 
  }
}
names(results)<-c("HR","lower95%CI","upper95%CI","p","feature") # HR=hazard ratio
print(results)

# Which var significant in the univariate analysis?
print("Variables which result to be significative in the univariate analysis: ")
print(results$feature[results$p<0.05])
sig_vars<-results$feature[results$p<0.05]

####################### DESCRIPTIVE ANALYSIS ################################################

# Steatoepatite
print("Steatoepatite")
summary(data$Steatoepatite)
summary(data.avg$Steatoepatite)

# HCV
print("HCV")
summary(data$HCV)
summary(data.avg$HCV)

# HBV
print("HBV")
summary(data$HBV)
summary(data.avg$HBV)

# potus
print("Potus")
summary(data$Potus.)
summary(data.avg$Potus.)

# albuminemia - many NA, little difference between pre-MI and post-MI
print("Albuminemia")
summary(data$Albuminemia)
boxplot(data$Albuminemia, main="Albuminemia")

summary(data.avg$Albuminemia)

boxplot(data.avg$Albuminemia, main="Albuminemia")

# creatinina - many NA, little difference between pre-MI and post-MI, boxplot gets higher
print("Creatinina")
summary(data$Creatininemia)
boxplot(data$Creatininemia, main="Creatinina")

summary(data.avg$Creatininemia)
boxplot(data.avg$Creatininemia, main="Creatinina")

# BINARY: almost everyone to 0
# inv_macro_port
print("inv_macro_port")
summary(data$inv_macro_port)
summary(data.avg$inv_macro_port)

# inv_macro_sovraep
print("inv_macro_sovraep")
summary(data$inv_macro_sovraep)
summary(data.avg$inv_macro_sovraep)

# EdmondsonGrading - almost everyone to 2, some to 3, 1 to 4
print("Edmondson grading")
summary(data$EdmondsonGrading)
summary(data.avg$EdmondsonGrading)

# MVI - many NA
print("MVI")
summary(data$MVI)
summary(data.avg$MVI)

# satellitosi - almost everyone to 0
print("satellitosi")
summary(data$satellitosi)
summary(data.avg$satellitosi)

# R - many NA, almost everyone to 0
print("R")
summary(data$R)
summary(data.avg$R)
