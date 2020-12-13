#!/usr/bin/env Rscript
######################

STUDY_DIR="/scratch/rjirsara/study-PDS"
PROJECT_DIR="/scratch/rjirsara/proj20-BrainAgeEval"
SCRIPTS_DIR="/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval"

######
### Prepare Data By Relabeling Identifiers
######

MASTER<-read.csv(paste0(SCRIPTS_DIR,"/GradientTreeBoost/feature-names.csv"),header=F)
CROSS<-t(MASTER) ; colnames(CROSS)<-t(MASTER)[1,] ; CROSS<-CROSS[-c(1),]
IDS<-as.data.frame(matrix(NA,nrow=0,ncol=2))
colnames(IDS)<-c("sub","ses")

######
### Prepare Data By Relabeling Identifiers
######

INDEX<-0
for (FILE in list.files(paste0(PROJECT_DIR,"/freesurf_extract/Data"), full.names=T)){
	INDEX=INDEX+1
	CONTENT<-read.csv(FILE,header=F)
	TEMP<-merge(MASTER,CONTENT, by="V1", sort = F)
	ROW<-as.numeric(t(TEMP)[2,]) ; CROSS<-as.data.frame(rbind(CROSS,ROW))
	CROSS[1:1118] <- lapply(CROSS[1:1118], as.numeric)
	SUB<-gsub("_ses","",unlist(strsplit(basename(FILE),"-"))[2])
	SES<-gsub(".csv","",unlist(strsplit(basename(FILE),"-"))[3])
	IDS[INDEX,]<-c(SUB,SES)
}

######
### Split Data By Gender & Save Spreadsheets
######

FULL<-cbind(IDS,CROSS) ; FULL$sub<-as.numeric(FULL$sub)
FINAL<-read.csv(paste0(PROJECT_DIR,"/Analysis/n834_DataFreeze_20201122.csv"))
FINAL<-FINAL[,c("sub","ses","sex")] ; FINAL<-merge(FINAL,FULL,by=c("sub","ses"))
FEMALES<-FINAL[which(FINAL$sex==2),] ; MALES<-FINAL[which(FINAL$sex==1),]
write.csv(FULL,paste0(PROJECT_DIR,"/Analysis/n",nrow(FULL),"_FS-Cross_Glasser_12122020.csv"),row.names=F)
write.csv(MALES,paste0(PROJECT_DIR,"/Analysis/n",nrow(MALES),"-males_FS-Cross_Glasser_12122020.csv"),row.names=F)
write.csv(FEMALES,paste0(PROJECT_DIR,"/Analysis/n",nrow(FEMALES),"-females_FS-Cross_Glasser_12122020.csv"),row.names=F)

######
### Prepare Data By Relabeling Identifiers
######

load("~/Desktop/brainageModels.RData")
FEMALES<-paste0(PROJECT_DIR,"/Analysis/n",nrow(FEMALES),"-females_FS-Cross_Glasser_12122020.csv")
MALES<-paste0(PROJECT_DIR,"/Analysis/n",nrow(MALES),"-males_FS-Cross_Glasser_12122020.csv")
FEMALES_IDS<-FEMALES[,1:2] ; FEMALES<-FEMALES[,-c(1:3)]
MALES_IDS<-MALES[,1:2] ; MALES<-MALES[,-c(1:3)]

######
### Run GTB Model and Save Brain Ages Within Master DataFreeze
######

FEMALES_IDS$brainage_GTB <- predict(mdl_agepred_female, as.matrix(FEMALES))
MALES_IDS$brainage_GTB <- predict(mdl_agepred_male, as.matrix(MALES))
GTB<-rbind(MALES_IDS,FEMALES_IDS)

FINAL<-read.csv(paste0(PROJECT_DIR,"/Analysis/n834_DataFreeze_20201122.csv"))
FINAL<-merge(FINAL,GTB,by=c("sub","ses"),all=T)
write.csv(FINAL,paste0(PROJECT_DIR,"/Analysis/n834_DataFreeze_20201212.csv"),row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
