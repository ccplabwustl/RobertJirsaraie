#!/usr/bin/env Rscript
######################

invisible(lapply(c("stringr","tidyr","plyr","dplyr","tidyverse","bnstruct","anticlust"), require, character.only=TRUE))
LABELS<-read.table("/scratch/jirsaraie/RobertJirsaraie/toolbox/bids_apps/atlases/atl-Gordon_Labels_Parcels.1D")
PROJ<-"/scratch/jirsaraie/proj22-NeuroMap"
ROOT<-"/scratch/jirsaraie/study-HCD"
TODAY<-format(Sys.time(), "%Y%m%d")

###### 
### Aggregate Quality Assurance Metrics
######

ANAT_QA<-NULL
ErinReid<-read.table(paste0(ROOT,"/audits/HCD_usable_list_6-29-22.1D"))
for (DIR in list.files(paste0(ROOT,"/apps/freesurfer"), full.names=T)){
	if (basename(DIR)=="fsaverage"){
		next
	}
	print(paste0("Euler Number: ", basename(DIR)))
	EULER<-read.table(paste0(DIR,"/stats/avg.orig.nofix.euler"),header=F)
	if (grepl(gsub("_MR","",basename(DIR)),ErinReid)){
		EReid<-"Include"
	} else {
		EReid<-"Exclude"
	}
	ANAT_QA<-rbind(ANAT_QA,as.data.frame(c(basename(DIR),EULER,EReid),col.names=c('sub','EulerNum','InclxSurf_ERied')))
}

FUNC_QA<-NULL ; SUB<-0
VARS<-c("sub","MisRegNodes","NumVolOrig","MeanDisOrig","MedianDisOrig","NumVolProc","MeanDisProc","MedianDisProc")
FUNC_QA<-setNames(data.frame(matrix(ncol=length(VARS), nrow=0)),VARS) 
for (DIR in list.files(paste0(ROOT,"/apps/xcp-fcon"), full.names=T)){
	print(paste0("Framewise Displacement: ", basename(DIR))) ; SUB<-SUB+1
	FUNC_QA[SUB,1]<-basename(DIR)
	if(!file.exists(paste0(DIR,"/problematic.1D"))){
    	FUNC_QA[SUB,2]<-NA
	} else {
		FUNC_QA[SUB,2]<-str_c(cbind(c(read.table(paste0(DIR,"/problematic.1D"),header=F)$V1)), collapse = ",")
	}
	FUNC_QA[SUB,3]<-length(read.table(paste0(ROOT,"/bids/",basename(DIR),"/func/FDfilt4.txt"))$V1)
	FUNC_QA[SUB,4]<-mean(read.table(paste0(ROOT,"/bids/",basename(DIR),"/func/FDfilt4.txt"))$V1)
	FUNC_QA[SUB,5]<-median(read.table(paste0(ROOT,"/bids/",basename(DIR),"/func/FDfilt4.txt"))$V1)
	FUNC_QA[SUB,6]<-length(read.table(paste0(ROOT,"/bids/",basename(DIR),"/func/FDfilt4_censored.txt"))$V1)
	FUNC_QA[SUB,7]<-mean(read.table(paste0(ROOT,"/bids/",basename(DIR),"/func/FDfilt4_censored.txt"))$V1)
	FUNC_QA[SUB,8]<-median(read.table(paste0(ROOT,"/bids/",basename(DIR),"/func/FDfilt4_censored.txt"))$V1)
}

QA<-merge(ANAT_QA,FUNC_QA,by=c("sub"),all=T) ; QA<-QA[,colSums(is.na(QA))<nrow(QA)] ; QA$sub<-gsub("_MR","",QA$sub)
for (SUB in QA$sub){
	if (file.exists(paste0(ROOT,"/bids/",SUB,"_MR/dwi/qc_mot_rel.1D"))){
		QA[which(QA$sub==SUB),"dwi_mot_rel"]<-as.numeric(read.table(paste0(ROOT,"/bids/",SUB,"_MR/dwi/qc_mot_rel.1D")))
		QA[which(QA$sub==SUB),"dwi_mot_abs"]<-as.numeric(read.table(paste0(ROOT,"/bids/",SUB,"_MR/dwi/qc_mot_abs.1D")))
	} else {
		QA[which(QA$sub==SUB),"dwi_mot_rel"]<-NA
		QA[which(QA$sub==SUB),"dwi_mot_abs"]<-NA	
	}
}

### Incl1MG_20220716: Select the First Decent Quality Timepoint
QA[which(QA$InclxSurf_ERied == "Exclude"),"Incl1MG_20220716"]<-"Exclude"
QA[which(QA$NumVolProc < 900),"Incl1MG_20220716"]<-"Exclude"
SUBIDS<-QA[is.na(QA$Incl1MG_20220716),"sub"]
for (SUBID in SUBIDS[grepl("_V1",SUBIDS)]){
	SUBSET<-QA[grep(unlist(strsplit(SUBID,"_"))[1],QA$sub),]
	INCL<-SUBSET[which(is.na(SUBSET[,"Incl1MG_20220716"])==TRUE)[1],"sub"]
	QA[which(QA$sub==INCL),"Incl1MG_20220716"]<-"Include"
}
QA[is.na(QA$Incl1MG_20220716),"Incl1MG_20220716"]<-"Exclude" 

### Incl2BL_20220816: Select the Heighest Quality Timepoint 
QA[which(QA$InclxSurf_ERied == "Exclude"),"Incl2BL_20220816"]<-"Exclude" 
QA$SUB<-gsub("HCD","",unlist(strsplit(gsub("_MR","",QA$sub),"_"))[grep("HCD",unlist(strsplit(gsub("_MR","",QA$sub),"_")))])
QA$SES<-gsub("V","",unlist(strsplit(gsub("_MR","",QA$sub),"_"))[grep("V",unlist(strsplit(gsub("_MR","",QA$sub),"_")))])
for (SUB in QA$SUB){ 
	SUBSET<-QA[which(QA$SUB==SUB),] 
	SUBSET<-SUBSET[which(SUBSET$InclxSurf_ERied=="Include"),]
	rSUBSET<-SUBSET[which(SUBSET$NumVolOrig > 1500),]
	INCL<-rSUBSET[which(rSUBSET$MeanDisOrig == min(rSUBSET[,"MeanDisOrig"])),"sub"]
	QA[which(QA$sub==INCL),"Incl2BL_20220816"]<-"Include"
}
QA[is.na(QA$Incl2BL_20220816),"Incl2BL_20220816"]<-"Exclude"

QA<-QA[,c(1,ncol(QA)-1,ncol(QA),2,4:11,3,12,13)] ; QA[,c(1:3)]<-lapply(QA[,c(1:3)],as.character) 
QA[,c(4:12)]<-lapply(QA[,c(4:12)],as.numeric) ; QA[,c(13:15)]<-lapply(QA[,c(13:15)],factor)

###### 
### Aggregate Demographics and Clinical Outcomes
######

RCParent<-read.csv(paste0(ROOT,"/audits/HCD_RedCap-Parent_2022-11-22.csv"))
RCParent$sub<-paste0(RCParent$subject, "_", RCParent$redcap_event)
RCParent<-RCParent[,c("sub","age","gender","internal_scale_pre","external_scale_pre","tot_prob_scale","internal_scale_sch","external_scale_sch","total_scale")]
RCChild<-read.csv(paste0(ROOT,"/audits/HCD_RedCap-Child_2022-11-22.csv"))
RCChild$sub<-paste0(RCChild$subject, "_", RCChild$redcap_event)
RCChild<-RCChild[,c("sub","age","gender","ysra_intnl_scale","ysra_extnl_scale","ysra_c_scale")]
RCTeen<-read.csv(paste0(ROOT,"/audits/HCD_RedCap-Teen_2022-11-11.csv"))
RCTeen$sub<-paste0(RCTeen$subject, "_", RCTeen$redcap_event)
RCTeen<-RCTeen[,c("sub","age","gender","asr_syn_intnl","asr_syn_extnl","asr_syn_sum")]
RC<-merge(RCParent,RCChild,by="sub",all=T) ; RC<-merge(RC,RCTeen,by="sub",all=T)
RC<-RC[grep("_V",RC$sub),] ; row.names(RC)<-NULL ; RC<-RC[which(RC$sub %in% QA$sub),]
names(RC)[4:6]<-c("inter_v1_parent","exter_v1_parent","pfact_v1_parent")
names(RC)[7:9]<-c("inter_v2_parent","exter_v2_parent","pfact_v2_parent")
names(RC)[12:14]<-c("inter_v2_child","exter_v2_child","pfact_v2_child")
names(RC)[17:19]<-c("inter_v3_child","exter_v3_child","pfact_v3_child")

MASTER<-merge(QA,RC,by="sub",all=T) 
MASTER$age<-rowMeans(MASTER[,grep("age",names(MASTER))],na.rm=T)
MASTER[which(MASTER$age == 0),"age"]<-NA ; MASTER$age.x<-NULL ; MASTER$age.y<-NULL
MASTER$gender<-rowMeans(MASTER[,grep("gender",names(MASTER))],na.rm=T) 
MASTER[which(MASTER$gender == 0),"gender"]<-NA ; MASTER$gender.x<-NULL ; MASTER$gender.y<-NULL
MASTER[which(MASTER$gender==2),"gender"]<-"Female"
MASTER[which(MASTER$gender==1),"gender"]<-"Male"
for (SUB in unique(MASTER$SUB)){
	GENDER<-MASTER[which(MASTER$SUB == SUB),"gender"][1]
	MASTER[which(MASTER$SUB == SUB),"gender"]<-GENDER
}

for (PATHOLOGY in c("inter","exter","pfact")){
	TEMP<-MASTER[,c(1,grep(PATHOLOGY,names(MASTER)))]
	TEMP1<-pivot_longer(TEMP, 2:5, names_to = "scale", values_to = "score")
	TEMP2<-TEMP1[!is.na(TEMP1$score),]
	for (SUB in unique(TEMP2$sub)){
		SUBDATA<-TEMP2[which(TEMP2$sub==SUB),]
		if (length(grep("v2",SUBDATA$scale))==2){
			TEMP2<-TEMP2[-which(TEMP2$sub==SUB & TEMP2$scale==paste0(PATHOLOGY,"_v2_parent")),]
		}
	}
	names(TEMP2)<-c("sub","pathology_scale",paste0(PATHOLOGY,"_score"))
	if (PATHOLOGY!="pfact"){
		TEMP2[,2]<-NULL
	} else {
		TEMP2<-TEMP2[,c(1,3,2)]
	}
	MASTER<-merge(MASTER, TEMP2, by=c("sub"), all=T)
}
MASTER$pathology_scale<-gsub("pfact_","",MASTER$pathology_scale)
MASTER[grep("_v",names(MASTER))]<-NULL
MASTER[,c("gender","pathology_scale")]<-lapply(MASTER[,c("gender","pathology_scale")],factor) 
MASTER[,c("age","inter_score","exter_score","pfact_score")]<-lapply(MASTER[,c("age","inter_score","exter_score","pfact_score")],as.numeric)

######
### Aggregate Anatomical Features
######

ANAT<-NULL
for (DIR in list.files(paste0(ROOT,"/apps/freesurfer"), full.names=T)){
	if (basename(DIR)=="fsaverage"){
		next
	}
	print(paste0("Structural Features: ", basename(DIR)))
	GORDON<-as.data.frame(t(read.csv(paste0(DIR,"/stats/gordon.csv"),header=F)))
	GORDON<-cbind(setNames(as.data.frame(c(basename(DIR))),"sub"),setNames(GORDON[2,],GORDON[1,]))
	names(GORDON)<-gsub("_parcel_","_parc_",names(GORDON))
	if(dim(GORDON)[2] == 2370){
		for (REMOVE in c("gauscurv","meancurv","foldindex","curvindex")){
			GORDON<-GORDON[,-c(grep(REMOVE, names(GORDON)))]
		}
	}
	if (is.null(ANAT)){
		ANAT<-GORDON
	} else {
		ANAT<-rbind(ANAT,GORDON)
	}
}
ANAT[,2:ncol(ANAT)]<-lapply(ANAT[,2:ncol(ANAT)],as.numeric) 
ANAT$sub<-gsub("_MR","",ANAT$sub) 
MASTER<-merge(MASTER,ANAT,by="sub")

###### 
### Aggregate Functional Features
######

FUNC_FCON<-NULL 
for (DIR in list.files(paste0(ROOT,"/apps/xcp-fcon"), full.names=T)[1173:1687]){
	print(paste0("FCON: ", basename(DIR)))
	FCON<-read.table(paste0(DIR,"/fcon_parcels.tsv"),header=F)
	FCON[,"HEADER"]<-paste0("edge_",FCON$V1,"x",FCON$V2,"_fcon")
	for (ROW in 1:nrow(FCON)){
		PARC1<-FCON[ROW,"V1"]
		LABEL1<-unlist(strsplit(LABELS[PARC1,],"x"))[2]
		PARC2<-FCON[ROW,"V2"]
		LABEL2<-unlist(strsplit(LABELS[PARC2,],"x"))[2]
		FCON[ROW,"LABELER"]<-paste0("edge_",LABEL1,"x",LABEL2,"_fcon")
	}
	NETWORKS<-NULL
	NODES<-unique(gsub(".*x", "", LABELS$V1))
	for (REGION1 in NODES){
		for (REGION2 in NODES[seq(which(NODES==REGION1),13)]){
			EDGE<-paste0("edge_",REGION1,'x',REGION2,"_fcon")
			STRENGTH<-mean(FCON[which(FCON$LABELER==EDGE),"V3"]) #ABSOLUTE FCON?
			NETWORKS<-rbind(NETWORKS,cbind(basename(DIR),EDGE,STRENGTH))
		}
	}
	NETWORKS<-pivot_wider(as.data.frame(NETWORKS),names_from=EDGE,values_from=STRENGTH)
	names(NETWORKS)[1]<-"sub" ; FUNC_FCON<-rbind(FUNC_FCON,NETWORKS)
	write.csv(FUNC_FCON,paste0(ROOT,"/zfunc-fcon_",dim(FUNC_FCON)[1],"x",dim(FUNC_FCON)[2],"_",TODAY,".csv"),row.names=F)
} 
FUNC_FCON$sub<-gsub("_MR","",FUNC_FCON$sub) 

FUNC_ALFF<-NULL
for (DIR in list.files(paste0(ROOT,"/apps/xcp-alff"), full.names=T)){
	print(paste0("ALFF: ", basename(DIR)))
	ALFF<-read.csv(paste0(DIR,"/alff_quantifyAtlas.csv"),header=T)
	HEADER<-str_replace(str_replace(names(ALFF),"glasser__Left","lh_L"),"glasser__Right","rh_R")
	names(ALFF)<-paste0(HEADER,"_ROI_alff") ; names(ALFF)[1]<-"sub" 
	FUNC_ALFF<-rbind.fill(FUNC_ALFF,ALFF)
}
names(FUNC_ALFF)<-gsub("_ROI","",gsub("mean_alff","roi",names(FUNC_ALFF)))
names(FUNC_ALFF)[1]<-"sub" ; FUNC_ALFF$sub<-gsub("_MR","",FUNC_ALFF$sub)

FUNC_REHO<-data.frame(matrix(ncol = ncol(FUNC_ALFF), nrow = 0)) 
for (DIR in list.files(paste0(ROOT,"/apps/xcp-reho"), full.names=T)){
	print(paste0("REHO: ", basename(DIR)))
	REHO<-read.csv(paste0(DIR,"/reho_quantifyAtlas.csv"),header=T)
	FUNC_REHO<-rbind(FUNC_REHO,REHO)
}
colnames(FUNC_REHO)<-paste0(gsub("mean_reho","roi",names(FUNC_REHO)),"_reho") 
names(FUNC_REHO)[1]<-"sub" ; FUNC_REHO$sub<-gsub("_MR","",FUNC_REHO$sub)

DUALREG<-read.table("/scratch/jirsaraie/RobertJirsaraie/toolbox/bids_apps/atlases/atl-Gordon_Labels_Parcels.1D")
FUNC_DUALREG<-data.frame(matrix(ncol=dim(DUALREG)[1]+1,nrow=0)) ; names(FUNC_DUALREG)<-DUALREG$V1
names(FUNC_DUALREG)[ncol(FUNC_DUALREG)]<-"sub" ; ROW<-0
for (DIR in list.files(paste0(ROOT,"/apps/xcp-dualreg"), full.names=T)){
	if (basename(DIR) == "atlas"){
		next
	}
	print(paste0("DUAL REG: ", basename(DIR))) ;  ROW<-ROW+1
	FUNC_DUALREG[ROW,"sub"]<-basename(DIR)
	for (FILE in list.files(paste0(DIR,"/stats"),full.names=T,pattern="cohesion")){
		NODE<-as.numeric(str_replace(strsplit(basename(FILE),"_")[[1]][1],"node-",""))
		FUNC_DUALREG[ROW,NODE]<-read.table(FILE)
	}
}
FUNC_DUALREG<-FUNC_DUALREG[,c(ncol(FUNC_DUALREG),1:333)] 
colnames(FUNC_DUALREG)<-paste0("roi_",names(FUNC_DUALREG),"_cohesion") 
names(FUNC_DUALREG)[1]<-"sub" ; FUNC_DUALREG$sub<-gsub("_MR","",FUNC_DUALREG$sub)

###### 
### Aggregate Diffusion Features
######

DWI<-as.data.frame(matrix(ncol=0,nrow=0))
for (FEAT in c("FA","MD","OD","ICVF")){
	INDEX<-0
	print(paste0("DTI-",FEAT)) 
	for (SUB in list.files(paste0(ROOT,"/apps/dti-fsl"),pattern="HCD")){
		INDEX<-INDEX+1 
		DWI[INDEX, "sub"]<-gsub("_MR","",SUB)
		if (FEAT == "FA" || FEAT == "MD"){
			DIR="dti-fsl"
		} else {
			DIR="dti-noddi"
		}
		for (PARC in 1:48){
			FILE<-paste0(ROOT,"/apps/",DIR,"/",SUB,"/FEAT-",FEAT,"_JHU-",PARC,".1D")
			if (file.exists(FILE)){
				VALUE<-as.numeric(read.table(paste0(FILE)))
			} else {
					VALUE<-NA
			}
			DWI[INDEX, paste0("roi_",PARC,"xJHU_",FEAT)]<-VALUE
		}
	}
}
DWI[,2:ncol(DWI)]<-lapply(DWI[,2:ncol(DWI)],as.numeric)

######
### Merge Demo, Clinical, QA, and Neuroimaging Data
######

MASTER<-merge(MASTER,ANAT,by=c("sub"),all=TRUE)
MASTER<-merge(MASTER,FUNC_FCON,by=c("sub"),all=TRUE) 
MASTER<-merge(MASTER,FUNC_ALFF,by=c("sub"),all=TRUE)
MASTER<-merge(MASTER,FUNC_REHO,by=c("sub"),all=TRUE)
MASTER<-merge(MASTER,FUNC_DUALREG,by=c("sub"),all=TRUE)
MASTER<-merge(MASTER,DWI,by=c("sub"),all=TRUE)
write.csv(MASTER,paste0(ROOT,"/data/freeze_",nrow(MASTER),"x",ncol(MASTER)-1,"_",TODAY,".csv"),row.names=F)

######
### Inclusion Criteria for the NeuroMap Project
######

MASTER<-read.csv(paste0(ROOT,"/data/freeze_1688x2339_20230524.csv"))
DATA<-MASTER[complete.cases(MASTER),] #Remove Missing
DATA<-DATA[which(DATA$age > 8),]  #Remove Younger than 8 years
DATA$EulerNum<-abs(DATA$EulerNum) ; INCL<-DATA[,c('sub','EulerNum','MeanDisOrig','dwi_mot_rel')] 
for (DIM in c('EulerNum','MeanDisOrig','dwi_mot_rel')){
	IQR<-IQR(INCL[,DIM]) ; MEDIAN<-median(INCL[,DIM])
	LABEL1=paste0(DIM,"_INCL1") 
	INCL[which(INCL[,DIM] < (MEDIAN+IQR)*1.0),LABEL1]<-1
	INCL[which(INCL[,DIM] >= (MEDIAN+IQR)*1.0),LABEL1]<-0
	LABEL2=paste0(DIM,"_INCL2") 
	INCL[which(INCL[,DIM] < ((MEDIAN+IQR)*1.5)),LABEL2]<-1
	INCL[which(INCL[,DIM] >= ((MEDIAN+IQR)*1.5)),LABEL2]<-0
	ggplot(INCL, aes(x=INCL[,DIM], colour=factor(INCL[,LABEL2]), fill=factor(INCL[,LABEL2]))) +
  	geom_histogram(bins=250,color="black") + theme(legend.position="none") +
  	scale_fill_manual(values = c("#D70040", "#228B22")) + xlab(DIM) +
  	ggtitle(paste0("1.5 * Interquartile Range - Threshold: ",round((MEDIAN+IQR)*1.5,digits=2))) 
  ggsave(paste0(PROJ,"/figures/Incl-",LABEL2,".png"),width=9)
}
INCL[,"INCLUDE1"]<-floor(rowMeans(INCL[,grep("INCL1",names(INCL))]))
INCL[,"INCLUDE2"]<-floor(rowMeans(INCL[,grep("INCL2",names(INCL))]))
INCLUDE<-INCL[,c("sub","INCLUDE1","INCLUDE2")]
DATA<-merge(DATA,INCLUDE,by="sub") 
DATA_MAIN<-DATA[which(DATA$INCLUDE2==1),]
for (SUB in unique(DATA_MAIN[,'SUB'])){
	if (length(which(DATA_MAIN[,'SUB']==SUB)) != 1){
			SUBSET<-DATA_MAIN[which(DATA_MAIN[,'SUB']==SUB),]
			for (REMOVE in SUBSET[order(SUBSET$EulerNum),"sub"][-1]){
				DATA_MAIN<-DATA_MAIN[-which(DATA_MAIN$sub == REMOVE),]
			}
	}
}
DATA_MAIN<-DATA_MAIN[,-c(grep("INCLUDE",names(DATA_MAIN)))]
DATA_SENS<-DATA[which(DATA$INCLUDE1==1),]
for (SUB in unique(DATA_SENS[,'SUB'])){
	if (length(which(DATA_SENS[,'SUB']==SUB)) != 1){
			SUBSET<-DATA_SENS[which(DATA_SENS[,'SUB']==SUB),]
			for (REMOVE in SUBSET[order(SUBSET$EulerNum),"sub"][-1]){
				DATA_SENS<-DATA_SENS[-which(DATA_SENS$sub == REMOVE),]
			}
	}
}
DATA_SENS<-DATA_SENS[,-c(grep("INCLUDE",names(DATA_SENS)))]

######
### Execute AntiClustering to Define Training and Test Splits
######

DATA_MAIN$anticlusters <- factor(anticlustering(
  DATA_MAIN[,c("age","inter_score","pfact_score","pfact_score")],
  categories = DATA_MAIN[,c("gender")],
  objective = "variance",
  K = 5
)) ; by(DATA_MAIN[,c("age","inter_score","pfact_score","pfact_score")], DATA_MAIN$anticlusters, function(x) round(colMeans(x),2))
DATA_MAIN<-DATA_MAIN[,c(1:21,ncol(DATA_MAIN),22:(ncol(DATA_MAIN)-1))]
write.csv(DATA_MAIN,paste0(PROJ,"/data/freeze_main_",nrow(DATA_MAIN),"x",ncol(DATA_MAIN)-1,"_",TODAY,".csv"),row.names=F)

DATA_SENS$anticlusters <- factor(anticlustering(
  DATA_SENS[,c("age","inter_score","pfact_score","pfact_score")],
  categories = DATA_SENS[,c("gender")],
  objective = "variance",
  K = 5
)) ; by(DATA_SENS[,c("age","inter_score","pfact_score","pfact_score")], DATA_SENS$anticlusters, function(x) round(colMeans(x),2))
DATA_SENS<-DATA_SENS[,c(1:21,ncol(DATA_SENS),22:(ncol(DATA_SENS)-1))]
write.csv(DATA_SENS,paste0(PROJ,"/data/freeze_sens_",nrow(DATA_SENS),"x",ncol(DATA_SENS)-1,"_",TODAY,".csv"),row.names=F)

######
### Create Feature Map Dictionary
######

DATA_MAIN<-read.csv(paste0(PROJ,"/data/freeze_main_956x2340_20230530.csv"),header=T) ; DATA_MAIN<-DATA_MAIN[,-c(1:22)] 
DF<-data.frame(matrix(NA,nrow=dim(DATA_MAIN)[2],ncol=0)) ; DF[,"features"]<-names(DATA_MAIN)

#Feature Type
DF[grep("_FA",DF$features),"feature_type"]<-"FA"
DF[grep("_MD",DF$features),"feature_type"]<-"MD"
DF[grep("_OD",DF$features),"feature_type"]<-"OD"
DF[grep("_ICVF",DF$features),"feature_type"]<-"ICVF"
DF[grep("_fcon",DF$features),"feature_type"]<-"FCON"
DF[grep("_alff",DF$features),"feature_type"]<-"ALFF"
DF[grep("_reho",DF$features),"feature_type"]<-"ReHo"
DF[grep("_cohesion",DF$features),"feature_type"]<-"NetCohesion"
DF[grep("_thickness_",DF$features),"feature_type"]<-"CortThick"
DF[grep("_volume_",DF$features),"feature_type"]<-"CortVolume"
DF[grep("_area_",DF$features),"feature_type"]<-"SurfArea"
DF[is.na(DF$feature_type),"feature_type"]<-"SubVolume"

#Network Labels
NET<-DF[grep("_cohesion",DF$features),"features"]
for (ROI in gsub("_cohesion","",gsub("roi_","",NET))){
  parc=unlist(strsplit(ROI,"x"))[1]
  network=unlist(strsplit(ROI,"x"))[2]
  DF[grep(paste0("_parc_",parc),DF$features),"network_label"]<-network
  DF[grep(paste0("roi_",parc,"x"),DF$features),"network_label"]<-network
  DF[grep(paste0(network,"x",network),DF$features),"network_label"]<-network
}
DF[grep(paste0("xJHU_"),DF$features),"network_label"]<-NA

#Color Codes
DF[grep("FA",DF$feature_type),"feature_color"]<-"#00ff00"
DF[grep("MD",DF$feature_type),"feature_color"]<-"#008000"
DF[grep("OD",DF$feature_type),"feature_color"]<-"#8fbc8f"
DF[grep("ICVF",DF$feature_type),"feature_color"]<-"#004c00"
DF[grep("ALFF",DF$feature_type),"feature_color"]<-"#ffbaba"
DF[grep("ReHo",DF$feature_type),"feature_color"]<-"#ff7b7b"
DF[grep("FCON",DF$feature_type),"feature_color"]<-"#a70000"
DF[grep("NetCohesion",DF$feature_type),"feature_color"]<-"#ff0000"
DF[grep("CortVolume",DF$feature_type),"feature_color"]<-"#58CCED"
DF[grep("CortThick",DF$feature_type),"feature_color"]<-"#072F5F"
DF[grep("SurfArea",DF$feature_type),"feature_color"]<-"#1261A0"
DF[grep("SubVolume",DF$feature_type),"feature_color"]<-"#3895D3"

#Modality Label
DF[which(DF$feature_type == "FA"),"modality"]<-"anat"
DF[which(DF$feature_type == "MD"),"modality"]<-"anat"
DF[which(DF$feature_type == "OD"),"modality"]<-"anat"
DF[which(DF$feature_type == "ICVF"),"modality"]<-"anat"
DF[which(DF$feature_type == "FCON"),"modality"]<-"func"
DF[which(DF$feature_type == "ALFF"),"modality"]<-"func"
DF[which(DF$feature_type == "ReHo"),"modality"]<-"func"
DF[which(DF$feature_type == "NetCohesion"),"modality"]<-"func"
DF[which(DF$feature_type == "CortThick"),"modality"]<-"dwi"
DF[which(DF$feature_type == "CortVolume"),"modality"]<-"dwi"
DF[which(DF$feature_type == "SurfArea"),"modality"]<-"dwi"
DF[which(DF$feature_type == "SubVolume"),"modality"]<-"dwi"

write.csv(DF,paste0(PROJ,"/data/featmap_",nrow(DF),"x",ncol(DF),"_",TODAY,".csv"),row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######