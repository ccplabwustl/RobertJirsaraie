#!/usr/bin/env Rscript
######################

library(gamm4) ; library(lme4) ; library(nlme) ; library(parameters) ; library(performance)
library(ggplot2) ; library(ggforce) ; library(corrplot) ; library(ggeffects) 
library(tidyr) ; library(plyr) ; library(dplyr) ; library(purrr) 
library(knitr) ; library(psych) ; library(broom.mixed) ; library(ggbeeswarm)
library(Rmisc) ; library(merTools) ; library(lm.beta)
library(lmerTest) ; library(gamm4) ; library(mgcv)

cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

normalize <- function(x) {
    return((x- min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

LmerBetas <- function(object) {
	sdy <- sd(getME(object,"y"))
	sdx <- apply(getME(object,"X"), 2, sd)
	sc <- fixef(object)*sdx/sdy
	se.fixef <- coef(summary(object))[,"Std. Error"]
	se <- se.fixef*sdx/sdy
	return(data.frame(stdcoef=sc, stdse=se))
}

######
### Read and Clean the Datasets
######

PATH<-"~/Box/Research/proj20-BrainAgeEval/audits"
PDS<-read.csv(paste0(PATH,"/n834_DataFreeze_20210414.csv"))
PDS<-PDS[,c(1:9,34:35,12,10:33,36:45)] ; names(PDS)[c(10,11)]<-c("hGTBgam","hGTBlmer") 
PDS[which(PDS$ses>3),"scanner"]<-"prisma" ; PDS[which(PDS$ses < 4),"scanner"]<-"trio" 
PDS[which(PDS$sex==1),"sex"]<-"males" ; PDS[which(PDS$sex==2),"sex"]<-"females"
PDS$scanner<-factor(PDS$scanner,levels=c("trio","prisma"))
names(PDS)[which(names(PDS)=="INR_FIRST")]<-"INR"
PDS$sex<-as.factor(PDS$sex)

HCP<-read.csv(paste0(PATH,"/n789_DataFreeze_20210414.csv"))
HCP[!complete.cases(HCP$Income),"INR"]<-NA ; HCP<-HCP[complete.cases(HCP$tDBN),] 
HCP[which(HCP$sex==1),"sex"]<-"males" ; HCP[which(HCP$sex==2),"sex"]<-"females"
HCP$sex<-as.factor(HCP$sex) 

######
### Append the 5-Cognitive Domains For Each Study
######

###PDS
cogPDS<-read.csv(paste0(PATH,"/PDS_Cognition.csv"))
names(cogPDS)[1]<-"sub" ; #cogPDS<-cogPDS[,c(1,grep("AgeAdjustedSS",names(cogPDS)))] #UnadjustedSS
COG<-pivot_longer(cogPDS, cols=2:ncol(cogPDS), names_to="domain", values_to="scores") 
COG$ses<-as.integer(substr(COG$domain,2,2)) ; COG<-COG[,c(1,4,2,3)]
COG$domain<-substr(gsub("_UnadjustedSS","_RAW",gsub("_AgeAdjustedSS","_NORM",COG$domain)),3,20)
COG<-pivot_wider(COG,names_from=domain,values_from=scores)
COG$sub<-as.numeric(gsub("L","",COG$sub))
NORMS<-grep("_NORM",names(COG))
for (SUB in unique(COG$sub)){
	COG[which(COG$sub==SUB),paste0(names(COG)[NORMS],"_TIV")]<-as.list(colMeans(COG[which(COG$sub==SUB),NORMS],na.rm=T))
}
COG<-as.data.frame(COG) ; COG<-COG[which(paste0(COG[,1],"x",COG[,2]) %in% paste0(PDS[,1],"x",PDS[,2])),]
PDS<-merge(PDS,COG,by=c("sub","ses"),all=T) ;
print("Create Time-Invarying Variables for Each Cogntive Domain")
for (SUB in unique(PDS$sub)){
	for (ROW in which(PDS$sub==SUB)){
		PDS[ROW,grep("_TIV",names(PDS))]<-colMeans(PDS[which(PDS$sub==SUB),grep("_TIV",names(PDS))],na.rm=T)
	}
}
PDS<-PDS[!is.na(rowMeans(PDS[,grep("_TIV",names(PDS))],na.rm=F)),]

#HCP
for (COL in c(17:18)){
	INDEX<-0
	for (FILE in list.files(paste0(PATH),recursive=T,full.names=T,pattern="HCD_Toolbox_Scored")){
		INDEX=INDEX+1
		if (INDEX == 1){
			cogHCP<-read.csv(FILE) 
		} else {
			COMP<-read.csv(FILE) ; cogHCP<-rbind(cogHCP,COMP)
		}
	}
	cogHCP<-cogHCP[which(cogHCP$visit=="V1"),] ; HCP$sub<-gsub("_V1_MR","",HCP$sub)
	cogHCP<-cogHCP[,c(1,4,COL)] ; names(cogHCP)<-c("sub","domain","score") ; cogHCP<-unique(cogHCP) #17 
	cogHCP<-as.data.frame(cogHCP %>% pivot_wider(names_from = domain, values_from = score))
	cogHCP<-cogHCP[,colSums(is.na(cogHCP))<nrow(cogHCP)] ; cogHCP$sub<-gsub("HCA","HCD",cogHCP$sub)
	cogHCP$Flanker<-rowMeans(cogHCP[,grep("Flanker",names(cogHCP))],na.rm=T)
	cogHCP$PictVocab<-rowMeans(cogHCP[,grep("Vocabulary",names(cogHCP))],na.rm=T)
	cogHCP$ListSort<-rowMeans(cogHCP[,grep("List Sorting",names(cogHCP))],na.rm=T)
	cogHCP$PictSeq<-rowMeans(cogHCP[,grep("Picture Sequence",names(cogHCP))],na.rm=T)
	cogHCP$PatComp<-rowMeans(cogHCP[,grep("Pattern Comparison",names(cogHCP))],na.rm=T)
	cogHCP$sub<-unlist(strsplit(cogHCP$sub,"_"))[grep("HCD",unlist(strsplit(cogHCP$sub,"_")))]
	cogHCP<-cogHCP[,c(1,31:35)] ; COGHCP<-cogHCP[which(cogHCP$sub %in% HCP$sub),] 
	if (COL == 17){
		names(COGHCP)[2:ncol(COGHCP)]<-paste0(names(COGHCP)[2:ncol(COGHCP)],"_RAW")
	} else if (COL == 18){
		names(COGHCP)[2:ncol(COGHCP)]<-paste0(names(COGHCP)[2:ncol(COGHCP)],"_NORM")
	}
	HCP<-merge(HCP,COGHCP,by=c("sub"),all=T)
}

PDS[c("INR_T1w","INR_FIRST.1","INR_LAST","INTLdim","EXTLdim","INTL_HBQ","EXTL_HBQ","TtlPrb_Z_sum")] <- list(NULL) 
PDS[c("TtlPrbHBQ_Z_sum","CDIC_TT","AvgMDDCore","AvgINTL","AvgEXTL","AvgINTL_HBQ","AvgEXTL_HBQ")] <- list(NULL) 
HCP[grep("brainage_tDBN",names(HCP))] <- list(NULL) ; HCP[c("Income","site","split")] <- list(NULL)
PDS[c("ManualQA","scanner")] <- list(NULL) ; PDS[grep("brainage",names(PDS))] <- list(NULL) 
PDS[which(PDS$ses > 3),"DATA"]<-"PDS_PRISMA" ; PDS[which(PDS$ses < 4),"DATA"]<-"PDS_TRIO"
HCP$DATA<-"HCP_TEST" ; PDS$DATA<-as.factor(PDS$DATA) ; HCP$DATA<-as.factor(HCP$DATA)
PDS[,paste0(names(PDS[,4:9]),"_ERROR")]<-PDS[,4:9]-PDS[,"age"] ; PRISMA<-PDS[which(PDS$DATA=="PDS_PRISMA"),]
HCP[,paste0(names(HCP[,3:8]),"_ERROR")]<-HCP[,3:8]-HCP[,"age"] ; TRIO<-PDS[which(PDS$DATA=="PDS_TRIO"),] 
PRISMA[,paste0(names(PRISMA)[32:37],"_NORM")]<-as.data.frame(lapply(PRISMA[,32:37], normalize))
TRIO[,paste0(names(TRIO)[32:37],"_NORM")]<-as.data.frame(lapply(TRIO[,32:37], normalize))
HCP[,paste0(names(HCP)[23:28],"_NORM")]<-as.data.frame(lapply(HCP[,23:28], normalize))
for (MOD in c("DBN_ERROR","rDBN_ERROR","tDBN_ERROR","GTB_ERROR","rGTB_ERROR","tGTB_ERROR")){
	M<-lmer(PRISMA[,which(names(PRISMA)==MOD)] ~ age + (1|sub), PRISMA)
	PRISMA[names(residuals(M)),gsub("ERROR","aERROR",MOD)]<-residuals(M)
	M<-lmer(TRIO[,which(names(TRIO)==MOD)] ~ age + (1|sub), TRIO)
	TRIO[names(residuals(M)),gsub("ERROR","aERROR",MOD)]<-residuals(M)
	M<-lm(HCP[,which(names(HCP)==MOD)] ~ age, HCP)
	HCP[names(M$residuals),gsub("ERROR","aERROR",MOD)]<-M$residuals
}
PRISMA[,paste0(names(PRISMA)[grep("aERROR",names(PRISMA))],"_NORM")]<-as.data.frame(lapply(PRISMA[,grep("aERROR",names(PRISMA))], normalize))
TRIO[,paste0(names(TRIO)[grep("aERROR",names(TRIO))],"_NORM")]<-as.data.frame(lapply(TRIO[,grep("aERROR",names(TRIO))], normalize))
HCP[,paste0(names(HCP)[grep("aERROR",names(HCP))],"_NORM")]<-as.data.frame(lapply(HCP[,grep("aERROR",names(HCP))], normalize))
COMBINE<-rbind.fill(TRIO,PRISMA,HCP) ; COMBINE$IDS<-paste0(COMBINE$sub,"x",COMBINE$ses)
for (ROW in 1:dim(COMBINE)[1]){
	COMBINE[ROW,"BRAINAGE_SD_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_ERROR",names(COMBINE))[1:6]]]),na.rm=T)
	COMBINE[ROW,"BRAINAGE_SD_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_ERROR_NORM",names(COMBINE))]]),na.rm=T)
	COMBINE[ROW,"DBN_SD_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_ERROR",names(COMBINE))[1:3]]]),na.rm=T)
	COMBINE[ROW,"DBN_SD_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_ERROR_NORM",names(COMBINE))][1:3]]),na.rm=T)
	COMBINE[ROW,"GTB_SD_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_ERROR",names(COMBINE))[4:6]]]),na.rm=T)
	COMBINE[ROW,"GTB_SD_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_ERROR_NORM",names(COMBINE))][4:6]]),na.rm=T)
	COMBINE[ROW,"BRAINAGE_aSD_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_aERROR",names(COMBINE))[1:6]]]),na.rm=T)
	COMBINE[ROW,"BRAINAGE_aSD_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_aERROR_NORM",names(COMBINE))]]),na.rm=T)
	COMBINE[ROW,"DBN_aSD_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_aERROR",names(COMBINE))[1:3]]]),na.rm=T)
	COMBINE[ROW,"DBN_aSD_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_aERROR_NORM",names(COMBINE))][1:3]]),na.rm=T)
	COMBINE[ROW,"GTB_aSD_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_aERROR",names(COMBINE))[4:6]]]),na.rm=T)
	COMBINE[ROW,"GTB_aSD_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_aERROR_NORM",names(COMBINE))][4:6]]),na.rm=T)
}
COMBINE$DATA<-factor(COMBINE$DATA,levels=c("PDS_TRIO","PDS_PRISMA","HCP_TEST"))
PRISMA<-COMBINE[which(COMBINE$DATA=="PDS_PRISMA"),] 
TRIO<-COMBINE[which(COMBINE$DATA=="PDS_TRIO"),] 
HCP<-COMBINE[which(COMBINE$DATA=="HCP_TEST"),]
QC<-read.csv(list.files(paste0(PATH),recursive=T,full.names=T,pattern="QC.csv"))
names(QC)[1]<-"sub" ; HCP<-merge(HCP,QC,by="sub")

######
### Aggregate Neuroimaging Files
######

###Neuroimaging Trajectories Across All Participants
NEURO_PATH<-"/Users/Jirsaraie/Box/Research/proj20-BrainAgeEval/audits"
FS_PDS<-read.csv(paste0(NEURO_PATH,"/n829_FS-Cross_20201231.csv"))
FS_PDS$sex<-NULL ; bPDS<-merge(PDS,FS_PDS,by=c("sub","ses")) 
bPDS<-bPDS[,-c(4:dim(PDS)[2])] ; bPDS$ses<-as.numeric(bPDS$ses) 
bPRISMA<-bPDS[which(bPDS$ses>3),] ; bPRISMA$ses<-NULL
bTRIO<-bPDS[which(bPDS$ses<4),] ; bTRIO$ses<-NULL
FS_HCP<-read.csv(paste0(NEURO_PATH,"/n789_FS-Cross_20210414.csv"))
FS_HCP$brainage_GTB<-NULL ; FS_HCP$sub<-gsub("_V1_MR","",FS_HCP$sub) 
bHCP<-merge(HCP,FS_HCP,by=c("sub")) ; bHCP$ses<-NULL ; HCP$ses<-NULL
bHCP<-bHCP[,-c(3:dim(HCP)[2])] 
bPRISMA$DATA<-"PDS_Prisma" ; bTRIO$DATA<-"PDS_Trio" ; bHCP$DATA<-"HCP_Test"
DATAFRAMES<-list(bTRIO,bPRISMA,bHCP) ; FINAL<-as.data.frame(matrix(NA,nrow=1,ncol=9))
FINAL[1,]<-c("feature","age_bin","DEV_RAW_SD","DEV_RAW_COV","DEV_Z_SD","DEV_Z_COV","DEV_MIN_SD","DEV_MIN_COV","data")
for (INDEX in 1:length(DATAFRAMES)){
	DATASET<-as.data.frame(DATAFRAMES[INDEX])
	DATASET[which(DATASET$age<=sort(DATASET$age)[round(nrow(DATASET)*0.25)]),"age_bin"]<-"YOUNG"
	DATASET[which(DATASET$age>=sort(DATASET$age)[round(nrow(DATASET)*0.75)]),"age_bin"]<-"OLD"
	DATASET[abs(scale(DATASET$age)) <= sort(abs(scale(DATASET$age)))[round(nrow(DATASET)*0.25)],"age_bin"]<-"MID"
	for (NEURO in 3:1120){
		for (AGE_BIN in c("YOUNG","MID","OLD")){
			DEV_RAW_SD<-sd(DATASET[which(DATASET$age_bin == AGE_BIN),NEURO],na.rm=TRUE)
			DEV_RAW_COV<-sd(DATASET[which(DATASET$age_bin == AGE_BIN),NEURO],na.rm=TRUE)/mean(DATASET[which(DATASET$age_bin == AGE_BIN),NEURO],na.rm=TRUE)	
			DATASET[,"TEMP_Z"]<-scale(DATASET[,NEURO])
			DEV_Z_SD<-sd(DATASET[which(DATASET$age_bin == AGE_BIN),"TEMP_Z"],na.rm=TRUE)
			DEV_Z_COV<-sd(DATASET[which(DATASET$age_bin == AGE_BIN),"TEMP_Z"],na.rm=TRUE)/mean(DATASET[which(DATASET$age_bin == AGE_BIN),"TEMP_Z"],na.rm=TRUE)
			DATASET[,"TEMP_MINMAX"]<-normalize(DATASET[,NEURO])
			DEV_MIN_SD<-sd(DATASET[which(DATASET$age_bin == AGE_BIN),"TEMP_MINMAX"],na.rm=TRUE)
			DEV_MIN_COV<-sd(DATASET[which(DATASET$age_bin == AGE_BIN),"TEMP_MINMAX"],na.rm=TRUE)/mean(DATASET[which(DATASET$age_bin == AGE_BIN),"TEMP_MINMAX"],na.rm=TRUE)
			FINAL<-rbind(FINAL,c(names(DATASET)[NEURO],AGE_BIN,DEV_RAW_SD,DEV_RAW_COV,DEV_Z_SD,DEV_Z_COV,DEV_MIN_SD,DEV_MIN_COV,DATASET$DATA[1]))
		}
	}
}
names(FINAL)<-FINAL[1,] ; FINAL<-FINAL[-1,]
FINAL[,3:8]<-lapply(FINAL[,3:8], as.numeric)
FINAL$age_bin<-factor(FINAL$age_bin,levels=c("YOUNG","MID","OLD"))
FINAL$data<-factor(FINAL$data,levels=c("PDS_Trio","PDS_Prisma","HCP_Test"))
FINAL[grep(".area",FINAL$feature),"FEATURE"]<-"AREA"
FINAL[grep("_volume",FINAL$feature),"FEATURE"]<-"VOLUME"
FINAL[grep("_thickness",FINAL$feature),"FEATURE"]<-"THICKNESS"
FINAL[which(is.na(FINAL$FEATURE)==TRUE),"FEATURE"]<-"SUBCORT"
PDSTrio<-FINAL[which(FINAL$data=="PDS_Trio"),]
PDSPrisma<-FINAL[which(FINAL$data=="PDS_Prisma"),]
HCPTest<-FINAL[which(FINAL$data=="HCP_Test"),]
rTRIO<-ddply(PDSTrio,"age_bin", summarize, MEAN_RAW_SD=mean(DEV_RAW_SD),MEAN_RAW_COV=mean(DEV_RAW_COV),MEAN_Z_SD=mean(DEV_Z_SD),MEAN_Z_COV=mean(DEV_Z_COV),MEAN_MIN_SD=mean(DEV_MIN_SD),MEAN_MIN_COV=mean(DEV_MIN_COV))
rPRISMA<-ddply(PDSPrisma,"age_bin", summarize, MEAN_RAW_SD=mean(DEV_RAW_SD),MEAN_RAW_COV=mean(DEV_RAW_COV),MEAN_Z_SD=mean(DEV_Z_SD),MEAN_Z_COV=mean(DEV_Z_COV),MEAN_MIN_SD=mean(DEV_MIN_SD),MEAN_MIN_COV=mean(DEV_MIN_COV))
rHCP<-ddply(HCPTest,"age_bin", summarize, MEAN_RAW_SD=mean(DEV_RAW_SD),MEAN_RAW_COV=mean(DEV_RAW_COV),MEAN_Z_SD=mean(DEV_Z_SD),MEAN_Z_COV=mean(DEV_Z_COV),MEAN_MIN_SD=mean(DEV_MIN_SD),MEAN_MIN_COV=mean(DEV_MIN_COV))
rTRIO$data<-"PDS_Trio" ; rPRISMA$data<-"PDS_Prisma" ; rHCP$data<-"HCP_Test" 
RESULTS<-rbind(rTRIO,rPRISMA,rHCP) 
FINAL<-merge(FINAL,RESULTS,by=c("age_bin","data"))

###Create Visuals of Reliability Differences in Age
bPDS<-bPDS[,-c(4:dim(PDS)[2])] ; bPDS$ses<-as.numeric(bPDS$ses) 
bTRIO<-bPDS[which(bPDS$ses<4),] ; bTRIO[,4:ncol(bTRIO)]<-scale(bTRIO[,4:ncol(bTRIO)])
bPRISMA<-bPDS[which(bPDS$ses>3),] ; bPRISMA[,4:ncol(bPRISMA)]<-scale(bPRISMA[,4:ncol(bPRISMA)])
LONG_PRISMA <-pivot_longer(bPRISMA,cols=4:ncol(bPRISMA),names_to="ROIs",values_to="Values")
LONG_TRIO <-pivot_longer(bTRIO,cols=4:ncol(bTRIO),names_to="ROIs",values_to="Values")
LONG_PRISMA$ses<-NULL ; LONG_PRISMA$DATA<-"PDS_Prisma"
LONG_TRIO$ses<-NULL ; LONG_TRIO$DATA<-"PDS_Trio"
FS_HCP$brainage_GTB<-NULL ; FS_HCP$sub<-gsub("_V1_MR","",FS_HCP$sub) 
bHCP<-merge(HCP,FS_HCP,by=c("sub")) ; bHCP$ses<-NULL
bHCP<-bHCP[,-c(3:dim(HCP)[2])] ; bHCP[,3:ncol(bHCP)]<-scale(bHCP[,3:ncol(bHCP)]) 
LONG_HCP <-pivot_longer(bHCP,cols=3:ncol(bHCP),names_to="ROIs",values_to="Values")
LONG_HCP$DATA<-"HCP_Test" ; LONG<-rbind(LONG_TRIO,LONG_PRISMA,LONG_HCP) 
LONG$DATA<-factor(LONG$DATA,levels=c("PDS_Trio","PDS_Prisma","HCP_Test"))
LONG[grep(".area",LONG$ROIs),"FEATURE"]<-"AREA"
LONG[grep("_volume",LONG$ROIs),"FEATURE"]<-"VOLUME"
LONG[grep("_thickness",LONG$ROIs),"FEATURE"]<-"THICKNESS"
LONG[which(is.na(LONG$FEATURE)==TRUE),"FEATURE"]<-"SUBCORT"

######
### Save Aggregated Datasets
######

write.csv(LONG,paste0(NEURO_PATH,"/NEURO_Traj.csv"),row.names=F)
write.csv(FINAL,paste0(NEURO_PATH,"/NEURO_Grp.csv"),row.names=F)
write.csv(COMBINE,paste0(NEURO_PATH,"/COMBINE.csv"),row.names=F)
write.csv(PRISMA,paste0(NEURO_PATH,"/PRISMA.csv"),row.names=F)
write.csv(TRIO,paste0(NEURO_PATH,"/TRIO.csv"),row.names=F)
write.csv(HCP,paste0(NEURO_PATH,"/HCP.csv"),row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######