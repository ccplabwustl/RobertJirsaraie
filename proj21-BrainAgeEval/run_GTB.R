#!/usr/bin/env Rscript
######################

STUDY_DIR="/scratch/rjirsara/study-PDS"
PROJECT_DIR="/scratch/rjirsara/proj20-BrainAgeEval"
SCRIPTS_DIR="/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval"

######
### Prepare Data By Relabeling Identifiers
######

MASTER<-read.csv(paste0(SCRIPTS_DIR,"/GradientTreeBoost/feature-names.csv"),header=F)
CROSS<-t(MASTER) ; colnames(CROSS)<-t(MASTER)[1,] ; CROSS<-CROSS[-c(1),] ; LONG<-CROSS

######
### Prepare Data By Relabeling Identifiers
######

FILES<-list.files(paste0(PROJECT_DIR,"/analysis/brainages_GTB"), full.names=T)
INDEX<-0 ; IDS<-as.data.frame(matrix(NA,nrow=0,ncol=2)) ; colnames(IDS)<-c("sub","ses")
for (FILE in FILES[lapply(FILES,function(x) length(grep("_long.csv",x,value=FALSE))) == 0]){
	INDEX=INDEX+1
	CONTENT<-read.csv(FILE,header=F)
	TEMP<-merge(MASTER,CONTENT, by="V1", sort = F)
	ROW<-as.numeric(t(TEMP)[2,]) ; CROSS<-as.data.frame(rbind(CROSS,ROW))
	CROSS[1:1118] <- lapply(CROSS[1:1118], as.numeric)
	SUB<-gsub("_ses","",unlist(strsplit(basename(FILE),"-"))[2])
	SES<-gsub(".csv","",unlist(strsplit(basename(FILE),"-"))[3])
	IDS[INDEX,]<-c(SUB,SES)
}
CROSS_FULL<-cbind(IDS,CROSS)

INDEX<-0 ; IDS<-as.data.frame(matrix(NA,nrow=0,ncol=2)) ; colnames(IDS)<-c("sub","ses")
for (FILE in list.files(paste0(PROJECT_DIR,"/analysis/brainages_GTB"), full.names=T,pattern="_long.csv")){
	INDEX=INDEX+1
	CONTENT<-read.csv(FILE,header=F)
	TEMP<-merge(MASTER,CONTENT, by="V1", sort = F)
	ROW<-as.numeric(t(TEMP)[2,]) ; LONG<-as.data.frame(rbind(LONG,ROW))
	LONG[1:1118] <- lapply(LONG[1:1118], as.numeric)
	SUB<-gsub("_ses","",unlist(strsplit(basename(FILE),"-"))[2])
	SES<-gsub("_long.csv","",unlist(strsplit(basename(FILE),"-"))[3])
	IDS[INDEX,]<-c(SUB,SES)
}
LONG_FULL<-cbind(IDS,LONG)

######
### Split Data By Gender & Save Spreadsheets
######

FINAL<-read.csv(paste0(PROJECT_DIR,"/analysis/n834_DataFreeze_20201213.csv"))
CROSS_FULL$sub<-as.numeric(CROSS_FULL$sub) ; FINAL_CROSS<-merge(FINAL,CROSS_FULL,by=c("sub","ses"))
LONG_FULL$sub<-as.numeric(LONG_FULL$sub) ; FINAL_LONG<-merge(FINAL,LONG_FULL,by=c("sub","ses")) 
MALES_CROSS<-FINAL_CROSS[which(FINAL_CROSS$sex==1),] ; MALES_LONG<-FINAL_LONG[which(FINAL_LONG$sex==1),]
FEMALES_CROSS<-FINAL_CROSS[which(FINAL_CROSS$sex==2),] ; FEMALES_LONG<-FINAL_LONG[which(FINAL_LONG$sex==2),]

######
### Save Datasets
######

write.csv(FINAL_LONG,paste0(PROJECT_DIR,"/analysis/n",nrow(LONG_FULL),"_FS-Long_20201231.csv"),row.names=F)
write.csv(FINAL_CROSS,paste0(PROJECT_DIR,"/analysis/n",nrow(CROSS_FULL),"_FS-Cross_20201231.csv"),row.names=F)
#write.csv(MALES_CROSS,paste0(PROJECT_DIR,"/analysis/n",nrow(MALES_CROSS),"-MALES_FS-Cross_20201231.csv"),row.names=F)
#write.csv(FEMALES_CROSS,paste0(PROJECT_DIR,"/analysis/n",nrow(FEMALES_CROSS),"-FEMALES_FS-Cross_20201231.csv"),row.names=F)
#write.csv(MALES_LONG,paste0(PROJECT_DIR,"/analysis/n",nrow(MALES_LONG),"-MALES_FS-Long_20201231.csv"),row.names=F)
#write.csv(FEMALES_LONG,paste0(PROJECT_DIR,"/analysis/n",nrow(FEMALES_LONG),"-FEMALES_FS-Long_20201231.csv"),row.names=F)

######
### Run GTB Model and Save Brain Ages Within Master DataFreeze
######

load("~/Desktop/brainageModels.RData")
for (FILE in list.files(paste0(PROJECT_DIR,"/analysis"),full.names=T)[!grepl("DataFreeze",list.files(paste0(PROJECT_DIR,"/analysis"),full.names=T))]){
	CONTENT<-read.csv(FILE) ; #print(names(CONTENT)[1:20]) ; print(FILE)
	FEMALES<-CONTENT[which(CONTENT$sex == 2),] ; IDS_FEMALES<-FEMALES[,c("sub","ses")] ; IDS_FEMALES$BRAINAGE<-0 
	FEMALES<-FEMALES[ , -which(names(FEMALES) %in% c("sub","ses","age","sex","SITE","SCANNER"))] ; colnames(FEMALES)<-NULL
	MALES<-CONTENT[which(CONTENT$sex == 1),] ; IDS_MALES<-MALES[,c("sub","ses")] ; IDS_MALES$BRAINAGE<-0 
	MALES<-MALES[ , -which(names(MALES) %in% c("sub","ses","age","sex","SITE","SCANNER"))] ; colnames(MALES)<-NULL
	IDS_FEMALES$BRAINAGE <- predict(mdl_agepred_female, as.matrix(FEMALES))
	IDS_MALES$BRAINAGE <- predict(mdl_agepred_male, as.matrix(MALES))
	GTB<-rbind(IDS_MALES,IDS_FEMALES) 
	names(GTB)[3]<-paste0("brainage_GTB_",gsub("-",'x',unlist(strsplit(basename(FILE),'_'))[2]))
	FINAL<-merge(FINAL,GTB,by=c("sub","ses"))
}

load("~/Desktop/ML-mdl-8-22y_toRobertDeanna_fullbrain_xgboost.RData")
for (FILE in list.files(paste0(PROJECT_DIR,"/analysis"),full.names=T)[!grepl("DataFreeze",list.files(paste0(PROJECT_DIR,"/analysis"),full.names=T))]){
	CONTENT<-read.csv(FILE) 
	FEMALES<-CONTENT[which(CONTENT$sex == 2),] ; IDS_FEMALES<-FEMALES[,c("sub","ses")] ; IDS_FEMALES$BRAINAGE<-0 
	FEMALES<-FEMALES[ , -which(names(FEMALES) %in% c("sub","ses","age","sex","SITE","SCANNER"))] ; colnames(FEMALES)<-NULL
	MALES<-CONTENT[which(CONTENT$sex == 1),] ; IDS_MALES<-MALES[,c("sub","ses")] ; IDS_MALES$BRAINAGE<-0 
	MALES<-MALES[ , -which(names(MALES) %in% c("sub","ses","age","sex","SITE","SCANNER"))] ; colnames(MALES)<-NULL
	IDS_FEMALES$BRAINAGE <- predict(mdl_agepred_female, as.matrix(FEMALES))
	IDS_MALES$BRAINAGE <- predict(mdl_agepred_male, as.matrix(MALES))
	GTB<-rbind(IDS_MALES,IDS_FEMALES) 
	names(GTB)[3]<-paste0("brainage_rGTB_",gsub("-",'x',unlist(strsplit(basename(FILE),'_'))[2]))
	FINAL<-merge(FINAL,GTB,by=c("sub","ses"))
}

#FINAL<-FINAL[ , -which(names(FINAL) %in% c("brainage_GTB"))]
write.csv(FINAL,paste0(PROJECT_DIR,"/analysis/n834_DataFreeze_20210101.csv"),row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
