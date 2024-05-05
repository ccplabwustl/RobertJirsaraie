#!/usr/bin/env Rscript
######################

invisible(lapply(c("stringr","tidyr","plyr","dplyr","tidyverse","corrplot","treemapify","lattice","tune","stringi","parameters","lme4","lmerTest","ggpubr"), require, character.only=TRUE))
ROOT<-"/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj22-NeuroMap/analysis"
DF<-read.csv(paste0(ROOT,"/20221215_UnweightXUnthresh/frozen-predict_20221215.csv"))
TODAY<-format(Sys.time(), "%Y%m%d") ; theme_set(theme_bw(base_size = 16))

######
### Prepare Data By Relabeling Identifiers
######

OUTCOMES<-read.csv("/scratch/jirsaraie/study-HCD/data/n789_DataFreeze_20210414.csv")
OUTCOMES<-OUTCOMES[,c("sub","age")] ; names(OUTCOMES)[1]<-"sub"
FUNC<-read.csv(paste0(ROOT,"/data/func-363x1092_20220404.csv"))
ANAT<-read.csv(paste0(ROOT,"/data/anat-365x1038_20220404.csv"))
for (SUB in FUNC$sub){
  SUBID<-unlist(strsplit(SUB,"_"))[1]
  FUNC[which(FUNC$sub == SUB),"age"]<-OUTCOMES[which(OUTCOMES$sub == SUBID),"age"][1]
}

DATA<-merge(ANAT,OUTCOMES,by="sub")
DATA[,2:1038]<-lapply(DATA[,2:1038],scale)
LONG<-pivot_longer(DATA,cols=2:1038,names_to="LABELS",values_to="VALUES") 

# Anatomical
LONG[grep("_area",LONG$LABELS),"FEAT_TYPE"]<-"AREA"
LONG[grep("_volume",LONG$LABELS),"FEAT_TYPE"]<-"VOLUME"
LONG[grep("_thickness",LONG$LABELS),"FEAT_TYPE"]<-"THICKNESS"
LONG[which(is.na(LONG$FEAT_TYPE)==TRUE),"FEAT_TYPE"]<-"SUBCORT"
LONG$FEAT_TYPE<-factor(LONG$FEAT_TYPE) ; LONG<-LONG[order(LONG$VALUES),]
ggplot(LONG, aes(x=age,y=VALUES, group=LABELS,color=FEAT_TYPE)) + 
  geom_smooth(method = "loess",se=F,size=0.2,alpha=0.75,fullrange=F) + 
  scale_color_manual(values=c("#0062ff","#f70505","#00991a","#ffa600")) + 
  theme_classic() 

# Functional
DATA<-merge(FUNC,OUTCOMES,by="sub")
DATA[,4:1093]<-lapply(DATA[,4:1093],scale)
LONG<-pivot_longer(DATA,cols=4:1093,names_to="LABELS",values_to="VALUES") 
LONG$FEAT_TYPE<-gsub(".*_","",LONG$LABELS) 
LONG$FEAT_TYPE<-factor(LONG$FEAT_TYPE) ; LONG<-LONG[order(LONG$VALUES),]
ggplot(LONG, aes(x=age,y=VALUES, group=LABELS,color=FEAT_TYPE)) + 
  geom_smooth(method = "loess",se=F,size=0.2,alpha=0.75,fullrange=F) + 
  scale_color_manual(values=c("#0062ff","#f70505","#00991a","#ffa600")) + 
  theme_classic() 

# Visualize the Inclusion SubSamples By Enrollment
MASTER<-read.csv("/scratch/jirsaraie/study-HCD/data/master_1645x19_20220901.csv")
MASTER<-MASTER[!is.na(MASTER$pathology_scale),]
MASTER$pathology_scale<-factor(MASTER$pathology_scale,levels=c("v1_parent","v2_parent","v2_child","v3_child"))
#MASTER<-MASTER[which(MASTER$Incl3RJ_20220821=="Include"),]
MASTER<-MASTER[order(MASTER$age),] ; MASTER$SUBJECT<- 0 
for (x in 1:length(unique(MASTER$SUB))){
  subid<-unique(MASTER$SUB)[x]
  MASTER[which(MASTER$SUB==subid),"SUBJECT"]<-x
}
MASTER[,c("SUBJECT","SES")]<-lapply(MASTER[,c("SUBJECT","SES")],factor)
MASTER<-MASTER[which(MASTER$InclxSurf_ERied=="Include"),]
ggplot() + 
  geom_line(MASTER,mapping=aes(x=age, y=SUBJECT, group=SUBJECT),size=0.5,alpha=0.5) +
  geom_point(MASTER,mapping=aes(x=age, y=SUBJECT,color=pathology_scale),size=1.5) + 
  scale_color_manual(values=c("#06029e","#00d9ff","#0400ff","#60029e")) + 
  theme_classic() 

###### 
### Aggregate Predictions
######

AggroData <- function(SUBLABEL){
  FILES<-list.files(ROOT, full.names=T, recursive=T, pattern="frozen-predict")
  PREDICT<-as.data.frame(matrix(nrow=0,ncol=0)) 
  PERFORM<-as.data.frame(matrix(nrow=0,ncol=3))
  INDEX<-0
  for (FILE in FILES[grep(SUBLABEL,FILES)]){
    INDEX<-INDEX+1
    TYPE<-unlist(strsplit(unlist(strsplit(FILE,"/"))[10],"-"))[2]
    DATA<-read.csv(FILE) ; names(DATA)[2]<-paste0("pfact_Pred-",TYPE)
    if (dim(PREDICT)[1] == 0){
      PREDICT<-DATA
    } else {
      PREDICT<-merge(PREDICT,DATA,by="sub")
    }
    MAE<-mean(read.csv(gsub("predict","perform",FILE))$mae)
    COD<-mean(read.csv(gsub("predict","perform",FILE))$cod)
    PERFORM[INDEX,]<-c(TYPE,MAE,COD)
  }
  names(PERFORM)<-c("LABEL","MAE","COD")
  return(list(STORE,PERFORM))
}

FEAT <- AggroData("20230206_FEAT")[[1]]
NET_PERFORM <- AggroData("20230206_NET")[[2]]
FEAT_PERFORM <- AggroData("20230206_FEAT")[[2]]

FOREST<-read.csv(paste0(ROOT,"/20221215_UnweightXUnthresh/SeeTheForest-Impurity.csv")) 
names(FOREST)[2]<-"LABEL" ; FOREST<-FOREST[,c("LABEL","feature_color")] 
FOREST<-FOREST[!duplicated(FOREST),]
FEAT_PERFORM<-merge(FEAT_PERFORM,FOREST,by="LABEL")

MULTI<-read.csv(paste0(ROOT,"/20221215_UnweightXUnthresh/frozen-perform_20230209.csv"))
MAE<-mean(MULTI[grepl("pfact",MULTI$X),"mae"])
COD<-mean(MULTI[grepl("pfact",MULTI$X),"cod"])
FEAT_PERFORM<-rbind(FEAT_PERFORM,c("MultiModal",MAE,COD,"#000000"))
NET_PERFORM<-rbind(NET_PERFORM,c("Global",MAE,COD))

###### 
### First Figure of Prediction Accuracy
######

FEAT_PERFORM$MAE<-as.numeric(FEAT_PERFORM$MAE)
FEAT_PERFORM$LABEL<-factor(FEAT_PERFORM$LABEL,levels=FEAT_PERFORM$LABEL[rev(order(FEAT_PERFORM$MAE))])
ggplot(data = FEAT_PERFORM) +
  geom_vline(xintercept = mean(FEAT_PERFORM$MAE), color = "#000000") +
  geom_label(aes(y=LABEL, x=MAE, label=round(MAE,digits=2),fill=LABEL,color=LABEL), size=6) +
  theme(legend.position="none")  +
  scale_fill_manual(values = c(FEAT_PERFORM$feature_color)) +
  scale_color_manual(values = c(rep("#ffffff", length(levels(FEAT_PERFORM$LABEL)))))
ggsave(paste0(ROOT,"/figures/F1A_FeatAccuracy.png"),width = 7)

NET_PERFORM$MAE<-as.numeric(NET_PERFORM$MAE)
NET_PERFORM$LABEL<-factor(NET_PERFORM$LABEL,levels=NET_PERFORM$LABEL[rev(order(NET_PERFORM$MAE))])
ggplot(data = NET_PERFORM) +
  geom_vline(xintercept = mean(NET_PERFORM$MAE), color = "#000000") +
  geom_label(aes(y=LABEL, x=MAE, label=round(MAE,digits=2),fill=MAE), size=6) +
  theme(legend.position="none")  +
  scale_fill_gradient(low="#c70000",high="#ffffff") 
ggsave(paste0(ROOT,"/figures/F1B_NetAccuracy.png"),width = 7)

###### 
### Second Figure of Prediction Bias
######

BIAS<-as.data.frame(matrix(nrow=0,ncol=2)) 
for (PRED in names(FEAT)[grep("Pred",names(FEAT))]){
  LABEL<-unlist(strsplit(PRED,"_"))[1]
  FEATURE<-unlist(strsplit(PRED,"-"))[2]
  SUBSET<-DF[,c("sub",paste0(LABEL,"_score"))]
  SUBSET<-merge(SUBSET,FEAT,by="sub")
  M<-lm(formula =  SUBSET[, paste0(LABEL, "_score")] ~ SUBSET[, PRED])
  COEF<-round(model_parameters(M,standardize = "refit")$Coefficient[2],digits=2)
  BIAS<-rbind(BIAS,c(FEATURE, COEF))
}
names(BIAS)<-c("LABEL","SLOPE")
FIX<-DF[,c(2,35,grep("pfact_predict",names(DF)))]
FIX[,"pfact_pred"]<-rowSums(FIX[,3:7],na.rm=T)
M<-lm(formula =  FIX[,"pfact_score"] ~ FIX[, "pfact_pred"])
COEF<-round(model_parameters(M,standardize = "refit")$Coefficient[2],digits=2)
BIAS<-rbind(BIAS,c("MultiModal", COEF))
BIAS<-merge(BIAS,FEAT_PERFORM[,c("LABEL","feature_color")],by="LABEL")

BIAS$SLOPE<-as.numeric(BIAS$SLOPE)
BIAS$LABEL<-factor(BIAS$LABEL,levels=BIAS$LABEL[order(BIAS$SLOPE)]) ; BIAS<-BIAS[order(BIAS$SLOPE),]
ggplot(data = BIAS) +
  geom_vline(xintercept = 1, color = "#000000") +
  geom_label(aes(y=LABEL, x=SLOPE, label=SLOPE,fill=LABEL,color=LABEL), size=6) +
  theme(legend.position="none")  +
  scale_fill_manual(values = c(BIAS$feature_color)) +
  scale_color_manual(values = c(rep("#ffffff", length(levels(FEAT_PERFORM$LABEL)))))
ggsave(paste0(ROOT,"/figures/F4S_PredBias.png"),width = 7)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######