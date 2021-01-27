#!/usr/bin/env Rscript
######################

library(stringr) ; library(dplyr)
TODAY<-format(Sys.time(), "%Y%m%d")
covaPath <- "/scratch/rjirsara/study-PDS/audits/n945_Predictors_Master_20201105.csv"
covaData<-read.csv(covaPath) ; names(covaData)[1]<-c("sub") ; names(covaData)[2]<-c("ses") ; covaData$T1w<-0
covaData[c(1:3)]<-lapply(covaData[c(1:3)],factor) ; covaData[,c(4:19)]<-lapply(covaData[,c(4:19)],as.numeric)
for (file in list.files("/scratch/rjirsara/proj20-BrainAgeEval/mass_N4correct",pattern="T1w.nii.gz")){
	subid<-gsub("sub-","",unlist(strsplit(file,"_",fixed=T))[1])
	sesid<-gsub("ses-","",unlist(strsplit(file,"_",fixed=T))[2])
	subid<-gsub("(?<![0-9])0+", "", subid, perl = TRUE)
	rownums<-which(covaData$sub == subid)
	rownum<-rownums[which(covaData[rownums,"ses"] == sesid)]
	covaData[rownum,19]<-1
	print(paste0("Sub: ",subid," Ses: ",sesid," Has Row Number: ",rownum))
}
covaData<-covaData[which(covaData$T1w==1),] ; covaData$T1w<-NULL

DBNPath <- "/scratch/rjirsara/proj20-BrainAgeEval/brainages_DBN/n834_DeepBrainNet.csv"
DBNData<-read.csv(DBNPath) ; data<-merge(covaData,DBNData,by=c("sub","ses"),all=T)

######
### Append Euler Number to Latest DataFreeze
######

data$EulerNumber<-NA
for (INDEX in 1:nrow(data)){
	SUB<-data[INDEX,"sub"] ; SES<-data[INDEX,"ses"] ; if (nchar(SUB) == 2){ SUB=paste0("0",SUB)}
	FILE<-paste0("/scratch/rjirsara/study-PDS/apps/freesurfer/sub-",SUB,"_ses-",SES,"/stats/avg.orig.nofix.euler")
	if (file.exists(FILE)){ 
		CONTENT<-read.table(FILE,header=F)
		data[INDEX,"EulerNumber"]<-as.numeric(CONTENT[1])
	}
}

write.csv(data,paste0("/scratch/rjirsara/proj20-BrainAgeEval/analysis/n",dim(data)[1],"_DataFreeze_",TODAY,".csv"),row.names=F)

######
### Append Manual QA Ratings of Raw T1w Scans
######

DATAFREEZE <- "/scratch/rjirsara/proj20-BrainAgeEval/analysis/n834_DataFreeze_20210101.csv"
DATA<-read.csv(DATAFREEZE) 
AUDIT <- "/scratch/rjirsara/study-PDS/audits/PDS.csv"
QA<-read.csv(AUDIT) ; QA$Subid_fMRI<-as.numeric(gsub("L","",QA$Subid_fMRI))

DATA$ManualQA<-NA
for (INDEX in 1:nrow(DATA)){
	SUB<-DATA[INDEX,"sub"]  ; SES<-DATA[INDEX,"ses"] 
	VALUE<-QA[which(QA$Subid_fMRI==SUB),SES+1]
	print(paste0(SUB,"x",SES,"x",INDEX,"x",VALUE))
	if (grepl(VALUE,"NO SCAN", fixed = TRUE)){
		print(paste0(SUB,"x",SES," Missing QA Rating"))
	
	} else if (is.na(unlist(strsplit(as.character(VALUE),"_"))[3]=="BAD")){
		print(paste0(SUB,"x",SES," QA Rating is Inclusionary"))
		DATA[INDEX,"ManualQA"]<-1
	} else {
		print(paste0(SUB,"x",SES," QA Rating is Exclusionary"))
		DATA[INDEX,"ManualQA"]<-0
	}
}
DATA$ManualQA<-as.factor(DATA$ManualQA)
write.csv(DATA,DATAFREEZE,row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
