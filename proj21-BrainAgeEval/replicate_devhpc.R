#!/usr/bin/env Rscript
######################

library(stringr) ; library(dplyr) ; library(tidyverse)
TODAY<-format(Sys.time(), "%Y%m%d")
PROJECT_DIR<-"/scratch/rjirsara/proj20-BrainAgeEval"
SCRIPTS_DIR<-"/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval"
covaPath <-paste0(PROJECT_DIR,"/replicate_devhcp/n789_RAW.csv")
data<-read.csv(covaPath) ; data[c(1,3:4)]<-lapply(data[c(1,3:4)],factor)

#DBNPath <- "/scratch/rjirsara/proj20-BrainAgeEval/brainages_DBN/n834_DeepBrainNet.csv"
#DBNData<-read.csv(DBNPath) ; data<-merge(covaData,DBNData,by=c("sub","ses"),all=T)

######
### Append Euler Number to Latest DataFreeze
######

data$EulerNumber<-NA
for (INDEX in 1:nrow(data)){
	print(paste0("WORKING ON INDEX:", INDEX))
	SUB<-data[INDEX,"sub"]
	FILE<-paste0(PROJECT_DIR,"/replicate_devhcp/EulerNumber/",SUB,".txt")
	if (file.exists(FILE)){ 
		CONTENT<-read.table(FILE,header=F)
		data[INDEX,"EulerNumber"]<-as.numeric(CONTENT[1])
	}
}

######
### Prepare Data By Relabeling Identifiers
######

MASTER<-read.csv(paste0(SCRIPTS_DIR,"/GradientTreeBoost/feature-names.csv"),header=F)
FILES<-list.files(paste0(PROJECT_DIR,"/replicate_devhcp/GlasserLabel"),full.names=T,pattern=".csv")
for (FILE in FILES){
	print(paste0("WORKING ON FILE: ", basename(FILE)))
	CONTENT<-read.csv(FILE,header=F)
	MASTER<-merge(MASTER,CONTENT, by="V1", sort = F,all=T)
	names(MASTER)[ncol(MASTER)]<-gsub(".csv","",basename(FILE))
	print(dim(MASTER)[1])
}
CROSS<-as.data.frame(t(MASTER)) ; colnames(CROSS)<-MASTER[,1]; CROSS<-CROSS[-c(1),]
CROSS[,1:1118]<-apply(CROSS[,1:1118], 2, function(x) as.numeric(as.character(x)))
CROSS$sub<-row.names(CROSS) ; DATA<-merge(data,CROSS,by=c("sub"),all=T)
write.csv(DATA,paste0(DIR_PROJECT,"/replicate_devhcp/n789_FS-Cross_20210207.csv"),row.names=F)

######
### Prepare Data By Relabeling Identifiers
######

load(paste0(SCRIPTS_DIR,"/GradientTreeBoost/brainageModels.RData"))
FEMALES<-DATA[which(DATA$sex==2),] ; IDS_FEMALES<-as.data.frame(FEMALES[,c("sub")]) ; IDS_FEMALES$BRAINAGE<-0 
FEMALES<-FEMALES[ , -which(names(FEMALES) %in% c("sub","age","sex","site","EulerNumber"))]
IDS_FEMALES$BRAINAGE <- predict(mdl_agepred_female, as.matrix(FEMALES))
MALES<-DATA[which(DATA$sex==1),] ; IDS_MALES<-as.data.frame(MALES[,c("sub")]) ; IDS_MALES$BRAINAGE<-0 
MALES<-MALES[ , -which(names(MALES) %in% c("sub","age","sex","site","EulerNumber"))]
IDS_MALES$BRAINAGE <- predict(mdl_agepred_male, as.matrix(MALES))

load(paste0(SCRIPTS_DIR,"/GradientTreeBoost/ML-mdl-8-22y_toRobertDeanna_fullbrain_xgboost.RData"))
IDS_FEMALES$RBRAINAGE<-0 ; IDS_FEMALES$RBRAINAGE <- predict(mdl_agepred_female, as.matrix(FEMALES)) 
names(IDS_FEMALES)<-c("sub","brainage_GTB","brainage_rGTB")
IDS_MALES$RBRAINAGE<-0 ; IDS_MALES$RBRAINAGE <- predict(mdl_agepred_male, as.matrix(MALES)) 
names(IDS_MALES)<-c("sub","brainage_GTB","brainage_rGTB")
AGGREGATE<-DATA[ ,which(names(DATA) %in% c("sub","age","sex","site","EulerNumber"))]
GTB<-rbind(IDS_MALES,IDS_FEMALES) ; AGGREGATE<-merge(AGGREGATE,GTB,by=c("sub"))
write.csv(DATA,paste0(DIR_PROJECT,"/replicate_devhcp/n789_DataFreeze_20210208.csv"),row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
