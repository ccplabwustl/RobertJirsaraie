#!/usr/bin/env Rscript
######################

lapply(c("ggplot2","tidyverse","dplyr","reshape","effectsize","effsize","lme4","devtools","ggeffects","gamm4","parameters","plyr","lmerTest"), require, character.only=TRUE)
TIME<-read.csv("https://raw.githubusercontent.com/ccplabwustl/RobertJirsaraie/master/proj21-DivergeReview/n19_MultimodTimeline_20220706.csv")
META<-read.csv("https://raw.githubusercontent.com/ccplabwustl/RobertJirsaraie/master/proj21-DivergeReview/n20_MultimodLit_20220701.csv")
META[which(META$modality=="Multimodal"),"modality"]<-"Multi" ; META$X<-NULL ; META$features_scaled<-NULL 
META$modality<-factor(gsub(" ","",META$modality)) ; META$domain<-factor(gsub(" ","",META$domain))

CONVERT_H2ETA <- function(H_STAT,DOF){
  F_STAT <- H_STAT/(3-1)
  ETA_SQAURED=(F_STAT * (3-1))/(F_STAT * (3-1) + (DOF-3))
  return(ETA_SQAURED)
}

normalize <- function(x) {
    return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

normalize_reverse <- function(x) {
    return((x-max(x,na.rm=T)) /(min(x,na.rm=T)-max(x,na.rm=T)))
}

META[which(META$domain=="cogntion"),"domain"]<-"cognition"
META[,c(6,7,10)]<-lapply(META[,c(6,7,10)],as.numeric) ; META[,c(1:5,8,9)]<-lapply(META[,c(1:5,8,9)],factor)

#####
### Create Figure of Brain Age Timeline
#####

TIME<-TIME[1:14,]
TIME$brainage<-TIME$brainage-TIME$mutlimodal
LONG<-pivot_longer(TIME,cols=c(2:3),names_to="CRITERIA",values_to="COUNT") 
ggplot(LONG,aes(year))+geom_bar(aes(weight=COUNT,fill=CRITERIA)) + 
  scale_color_manual(values = c("#3c84c7", "#355e3b")) +
  scale_fill_manual(values = c("#3c84c7", "#355e3b")) +
  scale_x_continuous(breaks = round(seq(min(LONG$year), max(LONG$year), by = 1),1)) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F1_BrainAgeTimeline.pdf",width=13,height=8)

#####
### Obtain Publicly Shared Data 
#####

### Richard2018 - https://osf.io/gwqmr/
COG<-read.table("/Users/Jirsaraie/Desktop/data/cogscores_anonymised.csv",header=T)
MODS<-read.table("/Users/Jirsaraie/Desktop/data/predAges_anonymised.csv",header=T)
DATA<-merge(MODS,COG,by="id") ; DATA$CP_sess<-NULL ; DATA$MMSE<-NULL
OUTPUT<-data.frame(matrix(ncol=2, nrow = 0))
CORRELATION<-list() ; INDEX<-0
for (MODEL in 4:14){
  INDEX<-INDEX+1
  DATA[,"TEMP_BAG"]<-DATA[,MODEL]-DATA[,"Age"] 
  DATA[,"TEMP_BAG"]<-residuals(lm(DATA[,"TEMP_BAG"]~DATA[,"Age"]))
  for (COGNITION in 38:91){
    CORR<-list(cor(DATA[,"TEMP_BAG"],DATA[,COGNITION],use = "complete.obs"))
    CORRELATION<-append(CORRELATION,CORR)
  }
  OUTPUT[INDEX,1]<-names(DATA)[MODEL]
  OUTPUT[INDEX,2]<-round(mean(abs(unlist(CORRELATION))),digits=2)
  CORRELATION<-NULL
}

### Dadi2020
GIT<-"https://raw.githubusercontent.com/RobertJirsaraie/empirical_proxy_measures/master/figure-3/inputs/imaging.csv"
DATA<-read.csv(GIT)
DATA$modality<-factor(DATA$modality)
for (MOD in levels(DATA$modality)){
  print(paste0(MOD,"x",mean(DATA[which(DATA$modality==MOD),"mae"])))
}

#####
### Convert the Effect Sizes Onto a Standardized Scale
#####

### Liem et al. 2017 contained 3 groups and a total of 1177 participants
### https://www.researchgate.net/post/Anyone-know-how-to-calculate-eta-squared-for-a-Kruskal-Wallis-analysis
for (ROW in which(META$study=="Liem2017")){
  H_STAT<-META[ROW,"utility_raw"]
  META[ROW,"utility_converted"]<-sqrt(CONVERT_H2ETA(H_STAT,724))
}

for (ROW in which(META$metric=="cohen-d")){
  META[ROW,"utility_converted"]<-d_to_r(META[ROW,"utility_raw"])
}

for (ROW in which(META$metric=="corr")){
  META[ROW,"utility_converted"]<-abs(META[ROW,"utility_raw"])
}

for (ROW in which(META$metric=="rsquared")){
  META[ROW,"utility_converted"]<-sqrt(META[ROW,"utility_raw"])
}

for (ROW in which(META$metric=="stdBetas")){
  META[ROW,"utility_converted"]<-abs(META[ROW,"utility_raw"])
}

META<-META[-which(META$study=="Engemann2020"),] ; row.names(DATA)<-NULL #TEMP
META[,c(6,7,10,11)]<-lapply(META[,c(6,7,10,11)],as.numeric)
META[,c(1:5,8,9)]<-lapply(META[,c(1:5,8,9)],factor)

#####
### Create Figure of Model Accuracy
#####

theme_set(theme_bw(base_size = 16))
CHECKPOINT<-META
META[which(META$modality=="ASL"),"modality"]<-"FUNC"
META$study<-factor(META$study,levels=unique(META$study))
MODALITY_ORDER<-c("Multi","T1w","T2w","DWI","FUNC","MEG")
META$modality<-factor(META$modality,levels=MODALITY_ORDER)
#META<-META[!is.na(META$modality),] ; META$modality<-factor(META$modality,levels=MODALITY_ORDER)
#theme(legend.position="none",aspect.ratio=1/1.25,panel.grid.minor=element_blank(),axis.title=element_blank())

#Remove Extra Rows For Studies that had more than 1 Developmental Outcome
for (REMOVE in 2:5){
  META<-META[-grep(paste0("_",REMOVE), META$studyxdomain),]
  row.names(META)<-NULL
  META$studyxdomain<-factor(META$studyxdomain)
}
tMETA<-rMETA
for (REMOVE in 2:5){
  tMETA<-tMETA[-grep(paste0("_",REMOVE), tMETA$studyxdomain),]
  row.names(tMETA)<-NULL
  tMETA$studyxdomain<-factor(tMETA$studyxdomain)
}

#Min/Max Scale the MAEs and Number of NeuroImaging Features
for (STUDY in unique(META$studyxdomain)){
  for (COL in c(6,7)){
    if (any(is.na(META[which(META$studyxdomain==STUDY),COL])) == FALSE){
      if (COL==6){
          META[which(META$studyxdomain==STUDY),paste0(names(META)[COL],"_scaled")]<-normalize(META[which(META$studyxdomain==STUDY),COL])  
        } else {
          META[which(META$studyxdomain==STUDY),paste0(names(META)[COL],"_scaled")]<-normalize_reverse(META[which(META$studyxdomain==STUDY),COL])  
      }
    }
  }
}

#Average Accruracy of all Models By Modality Within Each Study
F1<-NULL 
for (STUDY in unique(META$study)){
  SUBSET<-META[which(META$study==STUDY),]
  ADD<-ddply(SUBSET,"modality",summarize,accuracy_mae=mean(accuracy_mae),accuracy_mae_scaled=mean(accuracy_mae_scaled),study=STUDY) 
  F1<-rbind(F1,ADD)
}
ACCURACY<-ddply(F1,"modality",summarize,accuracy_mae=mean(accuracy_mae),accuracy_mae_scaled=mean(accuracy_mae_scaled),study="Aggregate") 

#Figure of Accuracy by Modality
F1$study<-factor(F1$study,levels=unique(F1$study))
F1$modality<-factor(F1$modality,levels=MODALITY_ORDER)
ACCURACY$modality<-factor(ACCURACY$modality,levels=MODALITY_ORDER)
ggplot() +
  geom_line(data=F1,aes(x=modality,y=accuracy_mae_scaled,group=study,color=study),size=0.5,alpha=0.5) + 
  geom_point(data=F1,aes(x=modality,y=accuracy_mae_scaled,group=study,color=study),size=2,alpha=0.75) +
  geom_line(data=ACCURACY,aes(x=modality,y=accuracy_mae_scaled,group=study),size=0.7,alpha=0.7,color="#000000") + 
  geom_point(data=ACCURACY,aes(x=modality,y=accuracy_mae_scaled),size=3,alpha=0.9,color="#000000") 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F1_Accuracy.pdf",width=8,height=7)

#Figure of Accuracy by Dimensonality
ggplot() +
  geom_smooth(data=META,aes(x=features_dim_scaled,y=accuracy_mae_scaled,group=study,color=study),method='lm',size=0.5,alpha=0.5,fullrange=F,se=F) + 
  geom_smooth(data=META,aes(x=features_dim_scaled,y=accuracy_mae_scaled),color="#000000",method='lm',size=2,alpha=0.75,fullrange=F,se=F) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F1B_Accuracy.pdf",width=8,height=7)

#Statistical Analyses
anova(lmer(accuracy_mae_scaled ~ modality + (1|study),data=META, REML=TRUE))
model_parameters(lmer(accuracy_mae_scaled ~ modality + (1|study),data=META, REML=TRUE), standardize="basic")
model_parameters(lmer(accuracy_mae_scaled ~ features_dim_scaled + (1|study),data=META, REML=TRUE), standardize="basic")

#####
### Create Figure of Model Utility
#####

rMETA<-CHECKPOINT[which(CHECKPOINT$Inclusion=="C"),]
rMETA[which(rMETA$modality=="ASL"),"modality"]<-"FUNC"
rMETA$study<-factor(rMETA$study,levels=unique(rMETA$study))
MODALITY_ORDER<-c("Multi","T1w","T2w","DWI","FUNC")
rMETA$modality<-factor(rMETA$modality,levels=MODALITY_ORDER)
rMETA<-rMETA[!is.na(rMETA$modality),] 
rMETA$modality<-factor(rMETA$modality,levels=MODALITY_ORDER)

#Min/Max Scale the MAEs and Number of NeuroImaging Features
for (STUDY in unique(rMETA$studyxdomain)){
  for (COL in c(6,7,11)){
    if (any(is.na(rMETA[which(rMETA$studyxdomain==STUDY),COL])) == FALSE){
      rMETA[which(rMETA$studyxdomain==STUDY),paste0(names(rMETA)[COL],"_scaled")]<-normalize(rMETA[which(rMETA$studyxdomain==STUDY),COL])  
    }
  }
}

#Average Utility of all Models by their Modality Within Each Study xdomain
F2<-NULL 
for (STUDY in unique(rMETA$study)){
  SUBSET<-rMETA[which(rMETA$study==STUDY),]
  LABEL<-rMETA[which(rMETA$study==STUDY)[1],"study"]
  ADD<-ddply(SUBSET,"modality",summarize,utility_converted=mean(utility_converted),utility_converted_scaled=mean(utility_converted_scaled),studyxdomain=STUDY,study=LABEL) 
  F2<-rbind(F2,ADD)
}
UTILITY<-ddply(F2,"modality",summarize,utility_converted=mean(utility_converted),utility_converted_scaled=mean(utility_converted_scaled),study="Aggregate") 

#Figure of Utility by Modality
F2$study<-factor(F2$study,levels=unique(F2$study))
F2$modality<-factor(F2$modality,levels=MODALITY_ORDER)
UTILITY$modality<-factor(UTILITY$modality,levels=MODALITY_ORDER)
ggplot() +
  geom_line(data=F2,aes(x=modality,y=utility_converted,group=studyxdomain,color=study),size=0.5,alpha=0.5) + 
  geom_point(data=F2,aes(x=modality,y=utility_converted,group=studyxdomain,color=study),size=2,alpha=0.75) +
  geom_line(data=UTILITY,aes(x=modality,y=utility_converted,group=study),size=0.7,alpha=0.7,color="#000000") + 
  geom_point(data=UTILITY,aes(x=modality,y=utility_converted),size=3,alpha=0.9,color="#000000") 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F2A_Utility.pdf",width=8,height=7)

#Figure of Utility by Dimensonality
ggplot() +
  geom_smooth(data=rMETA,aes(x=features_dim_scaled,y=utility_converted,group=study,color=study),method='lm',size=0.5,alpha=0.5,fullrange=F,se=F) + 
  geom_smooth(data=rMETA,aes(x=features_dim_scaled,y=utility_converted),color="#000000",method='lm',size=2,alpha=0.75,fullrange=F,se=F) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F2B_Utility.pdf",width=8,height=7)

#Statistical Analyses
anova(lmer(utility_converted ~ modality + (1|study/studyxdomain),data=rMETA, REML=TRUE))
model_parameters(lmer(utility_converted ~ modality + (1|study/studyxdomain),data=rMETA, REML=TRUE), standardize="basic")
model_parameters(lmer(utility_converted ~ features_dim_scaled + (1|studyxdomain),data=rMETA, REML=TRUE), standardize="basic")

#####
### Create Figure of Utility for Each Study
#####

rMETA$domain<-as.character(rMETA$domain)
rMETA[which(rMETA$domain=='AlcoholUse'),"domain"]<-"Alcohol Use"
rMETA[which(rMETA$domain=='alzheimers'),"domain"]<-"Alzheimer's Disease"
rMETA[which(rMETA$domain=='bipolar'),"domain"]<-"Bipolar Disorder"
rMETA[which(rMETA$domain=='cognition'),"domain"]<-"Cognition"
rMETA[which(rMETA$domain=='education'),"domain"]<-"Education"
rMETA[which(rMETA$domain=='fitness'),"domain"]<-"Physical Fitness"
rMETA[which(rMETA$domain=='mild_cog_impair'),"domain"]<-"Cognitive Impariments"
rMETA[which(rMETA$domain=='schizophrenia'),"domain"]<-"Schizophrenia"
rMETA[which(rMETA$domain=='sleepquality'),"domain"]<-"Sleep Quality"
rMETA[which(rMETA$domain=='subject_cog_impair'),"domain"]<-"Cognitive Complaints"
rMETA[which(rMETA$domain=='telomerlength'),"domain"]<-"Telomere Length"
rMETA[which(rMETA$domain=='BPSystolic'),"domain"]<-"Blood Pressure"
rMETA[which(rMETA$domain=='StrokeRisk'),"domain"]<-"Stroke Risk"
rMETA$domain<-factor(rMETA$domain,levels=c("Alcohol Use","Blood Pressure","Physical Fitness","Sleep Quality","Stroke Risk","Telomere Length","Alzheimer's Disease","Cognitive Complaints","Cognitive Impariments","Cognition","Education","Bipolar Disorder","Schizophrenia"))      

F3<-NULL 
for (STUDY in unique(rMETA$studyxdomain)){
  SUBSET<-rMETA[which(rMETA$studyxdomain==STUDY),]
  LABEL<-rMETA[which(rMETA$studyxdomain==STUDY)[1],"study"]
  DOMAIN<-rMETA[which(rMETA$studyxdomain==STUDY)[1],"domain"]
  ADD<-ddply(SUBSET,"modality",summarize,utility_converted=mean(utility_converted),utility_converted_scaled=mean(utility_converted_scaled),studyxdomain=STUDY,study=LABEL,domain=DOMAIN) 
  F3<-rbind(F3,ADD)
}

ggplot() +
  geom_line(data=F3,aes(x=modality,y=utility_converted,group=domain,color=domain),size=0.75,alpha=0.5) + 
  geom_point(data=F3,aes(x=modality,y=utility_converted,group=domain,color=domain),size=2,alpha=0.5) +
  facet_wrap(~study)
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F4A_UtilityByDomain.pdf",width=12,height=7) 

#####
### Create More Granular Figure of Model Utility Based on Developmental Outcomes
#####

#Average Utility of all Models by their Modality Within Each Study 
FSUPP<-NULL 
for (DOMAIN in unique(rMETA$domain)){
  SUBSET<-rMETA[which(rMETA$domain==DOMAIN),]
  ADD<-ddply(SUBSET,"modality",summarize,utility_converted=mean(utility_converted),utility_converted_scaled=mean(utility_converted_scaled),domain=DOMAIN) 
  FSUPP<-rbind(FSUPP,ADD)
}

FSUPP$domain<-factor(FSUPP$domain,levels=c("Alcohol Use","Blood Pressure","Physical Fitness","Sleep Quality","Stroke Risk","Telomere Length","Alzheimer's Disease","Cognitive Complaints","Cognitive Impariments","Cognition","Education","Bipolar Disorder","Schizophrenia"))      
fUTILITY<-ddply(FSUPP,"modality",summarize,utility_converted=mean(utility_converted),utility_converted_scaled=mean(utility_converted_scaled),domain="Aggregate") 
ggplot() +
  geom_line(data=FSUPP,aes(x=modality,y=utility_converted,group=domain,color=domain),size=0.5,alpha=0.5) + 
  geom_point(data=FSUPP,aes(x=modality,y=utility_converted,group=domain,color=domain),size=2,alpha=0.75) +
  geom_line(data=fUTILITY,aes(x=modality,y=utility_converted,group=domain),size=0.7,alpha=0.7,color="#000000") + 
  geom_point(data=fUTILITY,aes(x=modality,y=utility_converted),size=3,alpha=0.9,color="#000000") 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F4A_UtilityByDomain.pdf",width=8,height=7) 

ggplot() +
  geom_smooth(data=rMETA,aes(x=features_dim_scaled,y=utility_converted,group=domain,color=domain),method='lm',size=1.25,alpha=0.1,fullrange=F,se=F) + 
  geom_smooth(data=rMETA,aes(x=features_dim_scaled,y=utility_converted),color="#000000",method='lm',size=2.5,alpha=0.75,fullrange=F,se=F) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F4B_UtilityByDomain.pdf",width=8,height=7) 

#Statistical Analyses
TEMP<-rMETA
for (DOMAIN in unique(rMETA$domain)){
  TEMP<-rMETA[which(rMETA$domain == DOMAIN),]
  if (length(unique(TEMP$study)) == 1){ 
    print(paste0(""))
    #print(paste0(DOMAIN," - Model 1"))
    #print(anova(lm(utility_converted ~ modality,data=TEMP)))
    #print(model_parameters(lm(utility_converted ~ modality,data=TEMP), standardize="basic"))
    print(paste0(DOMAIN," - Model 2"))
    print(model_parameters(lm(utility_converted ~ features_dim_scaled,data=TEMP), standardize="basic"))
  } else {
    print(paste0(""))
    #print(paste0(DOMAIN," - Model 1"))
    #print(anova(lmer(utility_converted ~ modality+ (1|study),data=TEMP, REML=TRUE)))
    #print(model_parameters(lmer(utility_converted ~ modality + (1|study),data=TEMP, REML=TRUE), standardize="basic"))
    print(paste0(DOMAIN," - Model 2"))
    print(model_parameters(lmer(utility_converted ~ features_dim_scaled + (1|study),data=TEMP, REML=TRUE), standardize="basic"))
  }
}

#####
### Create Figure of the Relationship Between Accuracy and Utility
#####

ggplot() +
  geom_point(data=rMETA,aes(x=accuracy_mae_scaled,y=utility_converted,color=study),size=1.25,alpha=0.6) +
  geom_smooth(data=rMETA,aes(x=accuracy_mae_scaled,y=utility_converted),color="#000000",method='lm',size=2.5,alpha=0.75,se=F) + 
  geom_smooth(data=rMETA,aes(x=accuracy_mae_scaled,y=utility_converted,group=study,color=study),method='lm',size=1.25,alpha=0.6,se=F) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F3A_Assocation.pdf",width=9,height=7)

#Statistical Analyses
model_parameters(lmer(utility_converted ~ accuracy_mae_scaled + (1|study/studyxdomain),data=rMETA, REML=TRUE), standardize="basic")
TEMP<-NULL
for (STUDY in unique(rMETA$study)){
  TEMP<-rMETA[which(rMETA$study == STUDY),]
  if (length(unique(TEMP$studyxdomain)) == 1){ 
    print(paste0("")) ; print(paste0(STUDY))
    print(model_parameters(lm(utility_converted ~ accuracy_mae_scaled,data=TEMP), standardize="basic"))
  } else {
    print(paste0("")) ; print(paste0(STUDY))
    print(model_parameters(lmer(utility_converted ~ accuracy_mae_scaled + (1|studyxdomain),data=TEMP, REML=TRUE), standardize="basic"))
  }
}

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######

###F2: Scatterplots of Paradox By Dimesionality - Collective
LONG<-pivot_longer(META,cols=c(11,14),names_to="DOMAIN",values_to="POINTS") 
UTILITY<-ggpredict(lmer(utility_converted ~ features_dim_scaled + (1|study), data = META), terms = "features_dim_scaled")
ACCURACY<-ggpredict(lmer(accuracy_mae_scaled ~ features_dim_scaled + (1|study), data = META), terms = "features_dim_scaled")
STACK1<-as.data.frame(cbind("ACCURACY",ACCURACY$x,ACCURACY$predicted)) 
STACK2<-as.data.frame(cbind("UTILITY",UTILITY$x,UTILITY$predicted)) 
PARADOX<-NULL ; PARADOX<-rbind(STACK1,STACK2)
names(PARADOX)<-c("DOMAIN","X","Y") ; PARADOX[,c("X","Y")]<-lapply(PARADOX[,c("X","Y")], as.numeric)
ggplot() +
  geom_point(data=LONG,aes(x=features_dim_scaled,y=POINTS,group=DOMAIN,color=DOMAIN),size=2,alpha=0.7) + 
  geom_smooth(data=PARADOX,aes(x=X,y=Y,group=DOMAIN,color=DOMAIN),method='lm',size=3,alpha=0.7,fullrange=F,se=F) + 
  scale_fill_manual(values=c("#777B7B","#777B7B","#800000","#800000")) +
  scale_color_manual(values=c("#777B7B","#777B7B","#800000","#800000")) +
  facet_wrap(~domain_cluster) +
  theme(legend.position="none",aspect.ratio=1/1.25,panel.grid.minor=element_blank(),axis.title=element_blank()) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F2_Collective.pdf",width=7.4,height=6.8)

###F2: Scatterplots of Paradox By Dimesionality - Modality Specific
PARADOX<-NULL ; LONG<-NULL
for (MODE in c("anat", "dwi", "func")){
  TEMP<-META[which(META$modality == MODE),c(2,11,13,14)]
  TEMP1<-pivot_longer(TEMP,cols=c(2,4),names_to="DOMAIN",values_to="POINTS") 
  TEMP1[,"MODE"]<-MODE ; LONG<-rbind(LONG,TEMP1)
  UTILITY<-ggpredict(lmer(utility_converted ~ features_dim_scaled + (1|study), data = TEMP), terms = "features_dim_scaled")
  ACCURACY<-ggpredict(lmer(accuracy_mae_scaled ~ features_dim_scaled + (1|study), data = TEMP), terms = "features_dim_scaled")
  STACK1<-as.data.frame(cbind(MODE,"ACCURACY",ACCURACY$x,ACCURACY$predicted)) 
  STACK2<-as.data.frame(cbind(MODE,"UTILITY",UTILITY$x,UTILITY$predicted)) 
  PARADOX<-rbind(PARADOX,STACK1,STACK2)
}
names(PARADOX)<-c("MODE","DOMAIN","X","Y") ; PARADOX[,c("X","Y")]<-lapply(PARADOX[,c("X","Y")], as.numeric)
ggplot() +
  geom_point(data=LONG,aes(x=features_dim_scaled,y=POINTS,group=DOMAIN,color=DOMAIN),size=2,alpha=1) + 
  geom_smooth(data=PARADOX,aes(x=X,y=Y,group=DOMAIN,color=DOMAIN),method='lm',size=3,alpha=1,fullrange=F,se=F) + 
  scale_fill_manual(values=c("#777B7B","#777B7B","#A30000","#A30000")) +
  scale_color_manual(values=c("#777B7B","#777B7B","#A30000","#A30000")) +
  facet_wrap(~MODE,ncol=1) + 
  theme(legend.position="none",aspect.ratio=1/2.25,panel.grid.minor=element_blank(),axis.title=element_blank()) 
ggsave(paste0("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F2_ModalitySpecific.pdf"),width=7.4,height=6.8)

###F3: Histograms Regarding the Accuracy/Utility Paradox - Collective
GGPREDICT<-ggpredict(lmer(utility_converted ~ accuracy_mae_scaled + (1|studyxdomain), data = META), terms = "accuracy_mae_scaled")
ggplot() +
  geom_point(data=META,aes(x=accuracy_mae_scaled,y=utility_converted,group=study,color=study),size=2,alpha=0.7) + 
  geom_smooth(data=GGPREDICT,aes(x=x,y=predicted),color="#000000",method='lm',size=3,alpha=0.9,fullrange=F,se=F) + 
  scale_fill_manual(values=c("#d11141","#ffc425","#00b159","#00aedb","#f37735")) +
  scale_color_manual(values=c("#d11141","#ffc425","#00b159","#00aedb","#f37735")) +
  facet_wrap(~domain_cluster) +
  theme(legend.position="none",aspect.ratio=1/1.25,panel.grid.minor=element_blank(),axis.title=element_blank()) 
ggsave("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F3_Collective-scaledMAE.pdf",width=7.4,height=6.8)

###F3: Scatterplots of Paradox By Dimesionality - Modality Specific
PARADOX<-NULL ; LONG<-NULL
for (MODE in c("anat", "dwi", "func")){
  TEMP<-META[which(META$modality == MODE),c(1,2,7,11,13,14)]
  GGPREDICT<-ggpredict(lmer(utility_converted ~ accuracy_mae_scaled + (1|study), data = TEMP), terms = "accuracy_mae_scaled")
  STACK<-as.data.frame(cbind(MODE,GGPREDICT$x,GGPREDICT$predicted)) 
  PARADOX<-rbind(PARADOX,STACK)
}
names(PARADOX)<-c("MODE","X","Y") ; PARADOX[,c("X","Y")]<-lapply(PARADOX[,c("X","Y")], as.numeric)
sMETA<-META[which(META$modality==c("anat")),] 
fMETA<-META[which(META$modality==c("func")),] 
dMETA<-META[which(META$modality==c("dwi")),]
mMETA<-rbind(sMETA,fMETA,dMETA) ; names(mMETA)[which(names(mMETA)=="modality")]<-"MODE"
ggplot() +
  geom_point(data=mMETA,aes(x=accuracy_mae_scaled,y=utility_converted,group=MODE,color=study),size=2,alpha=1) + 
  geom_smooth(data=PARADOX,aes(x=X,y=Y),method='lm',color="#000000",size=3,alpha=1,fullrange=F,se=F) + 
  scale_fill_manual(values=c("#d11141","#ffc425","#00b159","#00aedb","#f37735")) +
  scale_color_manual(values=c("#d11141","#ffc425","#00b159","#00aedb","#f37735")) +
  facet_wrap(~MODE,ncol=1) + 
  theme(legend.position="none",aspect.ratio=1/2.25,panel.grid.minor=element_blank(),axis.title=element_blank()) 
ggsave(paste0("/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj21-DivergeReview/figures/F3_ModalitySpecific-scaledMAE.pdf"),width=7.4,height=6.8)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######