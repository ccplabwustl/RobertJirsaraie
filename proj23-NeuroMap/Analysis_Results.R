#!/usr/bin/env Rscript
######################

invisible(lapply(c("stringr","tidyr","plyr","dplyr","tidyverse","corrplot","lattice","stringi","parameters","lme4","lmerTest","ggplot2","effectsize"), require, character.only=TRUE))
ROOT<-"/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj22-NeuroMap/analysis/20230603_ADA"
FIGURES<-"/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj22-NeuroMap/figures/sensitivity"
TODAY<-format(Sys.time(), "%Y%m%d") 
theme_set(theme_bw(base_size = 16))
FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    if (is.null(x_scales) || is.null(x_scales)) {
        stop("X and Y scales required for facet_equal_wrap")
    }
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
  }
)
facet_wrap_equal <- function(...) {
  facet_super <- facet_wrap(...)
  ggproto(NULL, FacetEqualWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}

###### 
### Load & Prepare the Data
######

rm("fEXPLAIN") ; rm("nEXPLAIN") ; rm("fRESAMPL")
df = read.csv('/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj22-NeuroMap/data/freeze_main_956x2340_20230601.csv')
featmap = read.csv('/Users/Jirsaraie/Library/CloudStorage/Box-Box/Research/proj22-NeuroMap/data/featmap_2319x5_20230618.csv')
anat<-which(featmap$modality=="anat") 
dwi<-which(featmap$modality=="dwi") 
featmap[anat,"modality"]<-"dwi"
featmap[dwi,"modality"]<-"anat"
for (DIR in list.files(ROOT, full.names=T, recursive=F)){
  if (basename(DIR) == "Aggregate"){
    next
  }
  FILES<-list.files(DIR, full.names=T, recursive=T, pattern="csv")
  #EXPLAIN 
  if (basename(DIR) == "FEAT-Multimod"){
    EXPLAIN<-read.csv(FILES[grep("frozen-explain",FILES)])
    for (DIM in unique(EXPLAIN$dimension)){
      SUBSET<-as.data.frame(colMeans(EXPLAIN[which(EXPLAIN$dimension == DIM),-c(1:2)]))
      SUBSET[,"features"]<-row.names(SUBSET) ; row.names(SUBSET)<-NULL
      names(SUBSET)[1]<-"importance" ; SUBSET[,"dimensions"]<-DIM
      if (exists("fEXPLAIN")){
        fEXPLAIN<-rbind(fEXPLAIN,SUBSET)
      } else {
        fEXPLAIN<-SUBSET
      }
    }
  }
  if (basename(DIR) == "NET-Global"){
    EXPLAIN<-read.csv(FILES[grep("frozen-explain",FILES)])
    for (DIM in unique(EXPLAIN$dimension)){
      SUBSET<-as.data.frame(colMeans(EXPLAIN[which(EXPLAIN$dimension == DIM),-c(1:2)]))
      SUBSET[,"features"]<-row.names(SUBSET) ; row.names(SUBSET)<-NULL
      names(SUBSET)[1]<-"importance" ; SUBSET[,"dimensions"]<-DIM
      if (exists("nEXPLAIN")){
        nEXPLAIN<-rbind(nEXPLAIN,SUBSET)
      } else {
        nEXPLAIN<-SUBSET
      }
    }
  }
  #PERFORM
  SUBSET<-read.csv(FILES[grep("frozen-perform",FILES)])
  SUBSET[,"feature_type"]<-basename(DIR)
  if (exists("fPERFORM")){
    fPERFORM<-rbind(fPERFORM,SUBSET)
  } else {
    fPERFORM<-SUBSET
  }
  #PREDICT
  SUBSET<-read.csv(FILES[grep("frozen-predict",FILES)]) ; SUBSET[,"model"]<-basename(DIR)
  names(SUBSET)[1:4]<-c("sub","pfact_predict","inter_predict","exter_predict")
  if (exists("fPERDICT")){
    fPERDICT<-rbind(fPERDICT,SUBSET)
  } else {
    fPERDICT<-SUBSET
  }
  #RESAMPLE
  SUBSET<-read.csv(FILES[grep("frozen-resampl",FILES)])
  SUBSET[,"feature_type"]<-basename(DIR)
  if (exists("fRESAMPL")){
    fRESAMPL<-rbind(fRESAMPL,SUBSET)
  } else {
    fRESAMPL<-SUBSET
  }
}
featmap[grep("Retro",featmap$network_label),"network_label"]<-"NET-RetrosplTemporal"
fPERFORM[grep("Retro",fPERFORM$feature_type),"feature_type"]<-"NET-RetrosplTemporal"
fRESAMPL[grep("Retro",fRESAMPL$feature_type),"feature_type"]<-"NET-RetrosplTemporal"
fRESAMPL[grep("inter_score",fRESAMPL$dimension),"dimension"]<-"Internalized Symptoms"
fPERFORM[grep("inter_score",fPERFORM$dimension),"dimension"]<-"Internalized Symptoms"
fRESAMPL[grep("exter_score",fRESAMPL$dimension),"dimension"]<-"Externalized Symptoms"
fPERFORM[grep("exter_score",fPERFORM$dimension),"dimension"]<-"Externalized Symptoms"
fRESAMPL[grep("pfact_score",fRESAMPL$dimension),"dimension"]<-"General Psychopathology Factor"
fPERFORM[grep("pfact_score",fPERFORM$dimension),"dimension"]<-"General Psychopathology Factor"
fRESAMPL$dimension<-factor(fRESAMPL$dimension,levels=c("General Psychopathology Factor","Internalized Symptoms","Externalized Symptoms"))
FEAT_ORDER<-c("Multimod","CortThick","SurfArea","CortVolume","SubVolume","FCON","ALFF","ReHo","NetCohesion","FA","MD","ODI","ICVF")
suppressWarnings(dir.create(paste0(ROOT,"/Aggregate"), recursive = TRUE))
write.csv(nEXPLAIN,paste0(ROOT,"/Aggregate/nEXPLAIN.csv"),row.names=F)
write.csv(fEXPLAIN,paste0(ROOT,"/Aggregate/fEXPLAIN.csv"),row.names=F)
write.csv(fRESAMPL,paste0(ROOT,"/Aggregate/RESAMPLE.csv"),row.names=F)
write.csv(fPERFORM,paste0(ROOT,"/Aggregate/PERFORM.csv"),row.names=F)
write.csv(fPERDICT,paste0(ROOT,"/Aggregate/PREDICT.csv"),row.names=F)

#####
### Figure 1: Prediction Accuracy of Neural Properties
#####

FEAT<-fRESAMPL[grep("FEAT-", fRESAMPL$feature_type),] ; FEAT<-FEAT[order(FEAT$feature_type),]
FEAT[which(FEAT$feature_type=="FEAT-OD"),"feature_type"]<-"FEAT-ODI"
FEAT$feature_type<-factor(gsub("FEAT-","",FEAT$feature_type), levels=FEAT_ORDER) 

#Compute the Averages Based on Actual Performance
rm("FEAT_AVG")
for (DIM in unique(FEAT$dimension)){
  SUBSET1<-FEAT[grep(DIM, FEAT$dimension),]
  SUBSET1<-ddply(SUBSET1, "feature_type", summarize, MEAN_BTS=mean(mse), SD_BTS=sd(mse), dimension=DIM)
  SUBSET2<-fPERFORM[grep(DIM, fPERFORM$dimension),]
  SUBSET2<-SUBSET2[grep("FEAT-", SUBSET2$feature_type),]
  SUBSET2$feature_type<-factor(gsub("FEAT-","",SUBSET2$feature_type),levels=FEAT_ORDER)
  SUBSET2<-ddply(SUBSET2, "feature_type", summarize, MEAN_ACT=mean(mse), SD_ACT=sd(mse), dimension=DIM)
  SUBSET2<-SUBSET2[,-c(4)] ; SUBSET<-merge(SUBSET1,SUBSET2,by=c("feature_type"))
  if (exists("FEAT_AVG")){
    FEAT_AVG<-rbind(FEAT_AVG,SUBSET)
  } else {
    FEAT_AVG<-SUBSET
  }
}

#Add Color Codes
featmap[which(featmap$feature_type=="OD"),"feature_type"]<-"ODI"
for (TYPE in  unique(FEAT$feature_type)){
  if (TYPE != "Multimod"){
    FEAT[grep(TYPE, FEAT$feature_type),"feature_color"]<-featmap[grep(TYPE, featmap$feature_type)[1],"feature_color"]
    FEAT_AVG[grep(TYPE, FEAT_AVG$feature_type),"feature_color"]<-featmap[grep(TYPE, featmap$feature_type)[1],"feature_color"]
  } else {
    FEAT[grep(TYPE, FEAT$feature_type),"feature_color"]<-"#6F4E37"
    FEAT_AVG[grep(TYPE, FEAT_AVG$feature_type),"feature_color"]<-"#6F4E37"
  }
}

#Compute the Averages Based on Bootstrapped Samples
FEAT<-FEAT[which(FEAT$dimension=="General Psychopathology Factor"),]
for (GROUP in unique(FEAT$feature_type)){
  for (METRIC in 3:5){
    LABEL<-paste0(names(FEAT)[METRIC],"_avg")
    MEAN<-mean(FEAT[which(FEAT$feature_type==GROUP),METRIC])
    FEAT[which(FEAT$feature_type==GROUP),LABEL]<-MEAN
  }
}

#Plot Figure 1
FEAT<-FEAT[order(FEAT$feature_type),]
FIGURE1<-ggplot() + 
    geom_hline(yintercept = 0, color = "black", size = 0.25, alpha=0.5) + 
    geom_jitter(FEAT, width=0.35, size=0.5, alpha=0.5, mapping=aes(x=feature_type, y=cod, colour=feature_type)) + 
    geom_point(FEAT, mapping=aes(x=feature_type, y=cod_avg), size=7, shape=18, color="white") + 
    geom_point(FEAT, mapping=aes(x=feature_type, y=cod_avg), size=5, shape=18, color="black") + 
    scale_color_manual(values=unique(FEAT$feature_color)) +
    theme(strip.background = element_rect(colour = "black", fill = "black")) +
    theme(strip.text = element_text(color = "white", face="bold")) + 
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Neural Properties") + ylab("Model Accuracy (Coefficient of Determination)")
ggsave(plot=FIGURE1,filename=paste0(FIGURES,"/F1_",{TODAY},".png"),device="png",width=14,height=10.5,units='in')

#####
### Figure 2: Feature Importance of Multimodal Models
#####

#Prepare the Gini Impurity Data
F2<-merge(fEXPLAIN,featmap[,c("features","feature_type","feature_color")],by="features") 
F2<-F2[grep("pfact_score", F2$dimension),] ; F2$importance_percent<-(F2$importance)*100
F2<-merge(F2,ddply(F2, "feature_type", summarize, MEAN=mean(importance_percent)),by=c("feature_type"))
F2[which(F2$feature_type=="OD"),"feature_type"]<-"ODI"
F2$feature_type<-factor(F2$feature_type,levels=unique(F2[order(F2$MEAN,decreasing=TRUE),"feature_type"]))
F2<-F2[order(F2$MEAN,decreasing=TRUE),] 

#Plot Figure 2
FIGURE2<-ggplot() + 
    geom_hline(yintercept = 0, color = "black", size = 0.25, alpha=0.85) + 
    geom_jitter(F2, width=0.35, size=2.5,alpha=0.5, mapping=aes(x=feature_type, y=importance_percent, colour=feature_type)) + 
    geom_point(F2, mapping=aes(x=feature_type, y=MEAN), size=9, shape=18, color="white") + 
    geom_point(F2, mapping=aes(x=feature_type, y=MEAN), size=7, shape=18, color="black") + 
    scale_color_manual(values=unique(F2$feature_color)) +
    theme(legend.position="none") +
    xlab("Neural Properties") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Percent Contribution to Predictive Power (Gini Impurity * 100)")
ggsave(plot=FIGURE2,filename=paste0(FIGURES,"/F2_",{TODAY},".png"),device="png",width=14,height=10.5,units='in')
for (TYPE in unique(F2$feature_type)){
  SUBSET<-F2[which(F2$feature_type == TYPE),]
  FEAT<-SUBSET[order(SUBSET$importance,decreasing=T),"features"][1]
  PERCENT<-SUBSET[order(SUBSET$importance,decreasing=T),"importance_percent"][1]
  print(paste0(TYPE," most important feat is ",FEAT,"x",PERCENT))
}

#####
### Figure 3: Spatial Map of Important Neural Property
#####

#Plot Facet A of Figure 3
invisible(lapply(c("ggseg","ggseg3d","ggsegGordon","ggsegJHU"), require, character.only=TRUE))
F3<-merge(F2,featmap[,c("features","network_label","modality")],by="features")
F3[which(F3$network_label=="CinguloParietal"),"network_label"]<-"MedialParietal"
F3[which(F3$network_label=="NET-RetrosplTemporal"),"network_label"]<-"ParietoOccip"
GORDON<-as.data.frame(gordon$palette) 
GORDON$gordon<-row.names(GORDON) 
row.names(GORDON)<-NULL
for (ITERATION in GORDON$gordon){
  lNET<-unlist(strsplit(ITERATION,"_"))[1]
  cNET<-as.numeric(unlist(strsplit(ITERATION,"_"))[2])
  SUBSET<-F3[grep(lNET,F3$network_label),]
  SUBSET$features_number<-gsub("[^0-9.]", "", SUBSET$features)
  FEAT_NUM<-as.numeric(unique(SUBSET$features_number))[-1][cNET]
  KEEP_FCON<-SUBSET[which(SUBSET$feature_type=="FCON"),]
  KEEP<-SUBSET[which(SUBSET$features_number==FEAT_NUM),]
  KEEP<-KEEP[order(KEEP$importance,decreasing=TRUE),]
  GORDON[which(GORDON$gordon==ITERATION),"feature_color"]<-KEEP[1,"feature_color"]
  GORDON[which(GORDON$gordon==ITERATION),"feature_type"]<-c(KEEP[1,"feature_type"])
} ; gordon$data$type<-1
for (INDEX in 1:length(gordon$palette)){
  gordon$palette[INDEX]<-GORDON$feature_color[INDEX]
  gordon$data$type[INDEX]<-GORDON$feature_type[INDEX]
}
png(paste0(FIGURES,"/F3A_",TODAY,".png"),width=1008,height=288) 
plot(gordon) + labs("") + 
  theme(legend.position = "none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
dev.off() 

#Plot Facet B of Figure 3
F3B<-F3[which(F3$modality=="dwi"),] ; F3B$features_number<-as.numeric(gsub("[^0-9.]", "", F3B$features))
JHU<-as.data.frame(jhu$palette) 
for (INDEX in 1:dim(JHU)[1]){
  TRACT<-row.names(JHU)[INDEX]
  TRACT_INDEX<-which(jhu$data$region==TRACT)
  SUBSET<-F3B[which(F3B$features_number %in% TRACT_INDEX),]
  AVERAGE<-ddply(SUBSET,"feature_type",summarize,avg_value=mean(importance_percent))
  AVERAGE<-AVERAGE[order(AVERAGE$avg_value,decreasing=TRUE),]
  COLOR<-SUBSET[which(SUBSET$feature_type==AVERAGE$feature_type[1]),"feature_color"][1]
  TYPE<-SUBSET[which(SUBSET$feature_type==AVERAGE$feature_type[1]),"feature_type"][1]
  jhu$data[c(TRACT_INDEX),"color"]<-COLOR ; jhu$data[c(TRACT_INDEX),"type"]<-TYPE
} ; jhu$data$type<-factor(jhu$data$type)
FIGURE3B<-ggplot(data=jhu$data$geometry) + 
  geom_sf(aes(fill=jhu$data$type)) + 
  scale_fill_manual(values=c("#004c00","#8fbc8f","#008000","#00ff00")) + 
  theme(legend.position = "none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave(plot=FIGURE3B,filename=paste0(FIGURES,"/F3B_",TODAY,".png"),device="png",width=10.5,height=3,units='in')

#####
### Figure 4: Prediction Accuracy of Brain Networks
#####

#Compute the Averages Based on Actual Performance
rm("FEAT_AVG") ; rm("NET")
NET<-fRESAMPL[grep("NET-", fRESAMPL$feature_type),]
NET$feature_type<-gsub("NET-","", NET$feature_type) 
NET<-NET[order(NET$feature_type),]
for (DIM in unique(NET$dimension)){
  SUBSET1<-NET[grep(DIM, NET$dimension),]
  SUBSET1<-ddply(SUBSET1, "feature_type", summarize, MEAN_BTS=mean(cod), SD_BTS=sd(cod), dimension=DIM)
  SUBSET2<-fPERFORM[grep(DIM, fPERFORM$dimension),]
  SUBSET2<-SUBSET2[grep("NET-", SUBSET2$feature_type),] 
  SUBSET2$feature_type<-gsub("NET-","",SUBSET2$feature_type)
  SUBSET2<-ddply(SUBSET2, "feature_type", summarize, MEAN_ACT=mean(mse), SD_ACT=sd(mse), dimension=DIM)
  SUBSET2<-SUBSET2[,-c(4)] ; SUBSET<-merge(SUBSET1,SUBSET2,by=c("feature_type"))
  if (exists("FEAT_AVG")){
    FEAT_AVG<-rbind(FEAT_AVG,SUBSET)
  } else {
    FEAT_AVG<-SUBSET
  }
}

#Add Color Codes
FEAT_AVG$dimension<-factor(FEAT_AVG$dimension,levels=levels(NET$dimension))
PFACT<-FEAT_AVG[which(FEAT_AVG$dimension=="General Psychopathology Factor"),]
FEAT_AVG$feature_type<-factor(FEAT_AVG$feature_type,levels=c(PFACT[order(PFACT$MEAN_ACT),"feature_type"]))
NET$feature_type<-factor(NET$feature_type,levels=c(PFACT[order(PFACT$MEAN_ACT),"feature_type"]))
FEAT_AVG<-FEAT_AVG[order(FEAT_AVG$dimension),] ; FEAT_AVG<-FEAT_AVG[order(FEAT_AVG$feature_type),] 
colfunc<-colorRampPalette(c("#FFFC19","#FEF001","#FFCE03","#FD9A01","#FD6104","#FF2C05","#F00505","#6B1115"))
for (INDEX in seq(1,14)){
  NETWORK<-levels(FEAT_AVG$feature_type)[INDEX]
  FEAT_AVG[which(FEAT_AVG$feature_type==NETWORK),"feature_color"]<-colfunc(14)[INDEX]
}
PFACT_AVG<-FEAT_AVG[which(FEAT_AVG$dimension=="General Psychopathology Factor"),c("feature_type","MEAN_BTS","feature_color")]
NET<-NET[which(NET$dimension=="General Psychopathology Factor"),]
NET<-merge(NET,PFACT_AVG,by=c("feature_type"))
NET<-NET[order(NET$feature_type),]

#Plot Figure 4
FIGURE4<-ggplot() + 
    geom_hline(yintercept = 0, color = "black", size = 0.25, alpha=0.5) + 
    geom_jitter(NET, width=0.35, size=0.5, alpha=0.5, mapping=aes(x=feature_type, y=cod, colour=feature_type)) + 
    geom_point(NET, mapping=aes(x=feature_type, y=MEAN_BTS), size=7, shape=18, color="white") + 
    geom_point(NET, mapping=aes(x=feature_type, y=MEAN_BTS), size=5, shape=18, color="black") + 
    scale_color_manual(values=unique(NET$feature_color)) +
    theme(strip.background = element_rect(colour = "black", fill = "black")) +
    theme(strip.text = element_text(color = "white", face="bold")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="none") +
    xlab("Brain Networks") + ylab("Model Accuracy (Coefficient of Determination)")
ggsave(plot=FIGURE4,filename=paste0(FIGURES,"/F4_",{TODAY},".png"),device="png",width=14,height=10.5,units='in')

#####
### Figure 5: Heat map of Important Networks
#####

#Aggregate and Prepare the Gini Impurities
names(FEAT_AVG)[1]<-"network_label" 
Rfeatmap<-featmap[!is.na(featmap$network_label),]
Rfeatmap$network_label<-gsub("NET-","",Rfeatmap$network_label)
pfact_EXPLAIN<-nEXPLAIN[which(nEXPLAIN$dimensions=="exter_score"),] 
F5<-merge(pfact_EXPLAIN,Rfeatmap,by="features") 
F5<-F5[,c("features","network_label","feature_type","importance")]
rNET<-NET[,c("feature_type","feature_color")] 
names(rNET)[1]<-"network_label"
rNET<-rNET[!duplicated(rNET),]
F5$importance<-(F5$importance*100)
F5<-merge(F5,rNET,by=c("network_label")) 

#Compute the Averages Across Groups of Features
FINAL=as.data.frame(matrix(NA,nrow=0,ncol=4))
names(FINAL)<-c("NETWORK","TYPE","IMPURITY","COLOR")
for (TYPE in levels(factor(F5$feature_type))){
  SUBSET<-F5[which(F5$feature_type==TYPE),]
  for (GROUP in unique(SUBSET$network_label)){
    SUM<-sum(SUBSET[which(SUBSET$network_label==GROUP),"importance"])
    COLOR<- SUBSET[which(SUBSET$network_label==GROUP),"feature_color"][1]
    ADD_ROW<-as.data.frame(t(c(GROUP,TYPE,SUM,COLOR)))
    names(ADD_ROW)<-names(FINAL) 
    FINAL<-rbind(FINAL,ADD_ROW)
  }
} 
FINAL$IMPURITY<-round(as.numeric(FINAL$IMPURITY),digits=0)
FINAL$NETWORK<-factor(FINAL$NETWORK,levels=rev(unique(NET$feature_type)))
FINAL$TYPE<-factor(FINAL$TYPE,levels=c("CortThick","CortVolume","ReHo","NetCohesion","SurfArea","ALFF","FCON"))

#Plot Figure 5
theme_set(theme_classic(base_size = 16))
FIGURE5<-ggplot(data=FINAL,aes(x=TYPE,y=NETWORK,fill=NETWORK))+
  geom_tile(aes(alpha=IMPURITY),color="black")+
  geom_text(aes(label = sprintf("%0.0f",IMPURITY)),size = 5)+
  scale_alpha(range = c(0.25, 1))+
  xlab("Neural Properties") + ylab("Brain Networks") + 
  theme(strip.background = element_rect(colour = "black", fill = "black")) +
  theme(strip.text = element_text(color = "white", face="bold")) + 
  scale_fill_manual(values=rev(colfunc(13))) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), legend.position = "none") 
ggsave(plot=FIGURE5,filename=paste0(FIGURES,"/F5_",TODAY,".png"),device="png",width=14,height=10.5,units='in')

#####
### Figure 6: Spatial Map of Key Brain Networks 
#####

#Aggrogate the Relevant Packages
invisible(lapply(c("ggseg","ggseg3d","ggsegGordon","ggsegJHU"), require, character.only=TRUE))
Rfeatmap<-featmap[!is.na(featmap$network_label),]
Rfeatmap$network_label<-gsub("NET-","",Rfeatmap$network_label)
F6<-merge(nEXPLAIN,Rfeatmap[,c("features","feature_type","network_label")],by="features") 
F6<-F6[-which(F6$feature_type=="FCON"),] 
F6$importance<-F6$importance*100
F6$network_label<-as.character(F6$network_label)
F6[which(F6$network_label=="CinguloParietal"),"network_label"]<-"MedialParietal"
F6[which(F6$network_label=="RetrosplTemporal"),"network_label"]<-"ParietoOccip"
GORDON<-as.data.frame(gordon$palette) 
GORDON$gordon<-row.names(GORDON) 
row.names(GORDON)<-NULL

#Compute the Average Impurities Per Group
for (DIM in unique(F6$dimensions)){
  for (FEAT in unique(F6$feature_type)){
    for (NET in unique(F6$network_label)){
      SUBSET1<-F6[which(F6$dimensions==DIM),]
      SUBSET2<-SUBSET1[which(SUBSET1$feature_type==FEAT),]
      SUBSET3<-SUBSET2[which(SUBSET2$network_label==NET),]
      WEIGHT<-round(sum(SUBSET3[,"importance"]),digits=2)
      LABEL<-gsub("score",FEAT,DIM)
      GORDON[grep(NET,GORDON$gordon),gsub("score",FEAT,DIM)]<-WEIGHT
    }
  }
}
GORDON<-pivot_longer(GORDON, 3:20, names_to = "labels", values_to = "weights")
GORDON$feature_type <- sapply(strsplit(GORDON$labels, "_"), function(x) x[2])
GORDON$dimensions <- sapply(strsplit(GORDON$labels, "_"), function(x) x[1])
GORDON<-GORDON[,-c(1)] ; GORDON$labels<-NULL ; names(GORDON)[1]<-"region"
for (FEAT in unique(F6$feature_type)){
  for (DIM in unique(F6$dimension)){
    GORDON<-rbind(GORDON,c("NEUTRAL",NA,FEAT,gsub("_score","",DIM)))
  }
} ; GORDON$weights<-as.numeric(GORDON$weights)

#Perfect the Figure Lables
GORDON<-as.data.frame(GORDON)
GORDON[which(GORDON$feature_type=="CortThick"),"feature_type"]<-"Cortical Thickness"
GORDON[which(GORDON$feature_type=="CortVolume"),"feature_type"]<-"Cortical Volume"
GORDON[which(GORDON$feature_type=="SurfArea"),"feature_type"]<-"Surface Area"
GORDON[which(GORDON$feature_type=="NetCohesion"),"feature_type"]<-"Network Cohesion"
GORDON$feature_type<-factor(GORDON$feature_type,levels=c("Cortical Thickness","Cortical Volume","ReHo","Network Cohesion","Surface Area","ALFF"))
GORDON$dimensions<-as.character(GORDON$dimensions)
GORDON[which(GORDON$dimensions=="inter"),"dimensions"]<-"internalizing"
GORDON[which(GORDON$dimensions=="exter"),"dimensions"]<-"externalizing"
GORDON[which(GORDON$dimensions=="pfact"),"dimensions"]<-"p-factor"
GORDON$dimensions<-factor(GORDON$dimensions,levels=c("p-factor","internalizing","externalizing"))
GORDON<-as.tibble(GORDON)

#Plot Figure 6
FIGURE6 <- GORDON %>%
  group_by(feature_type) %>%
    ggplot() +
      geom_brain(atlas = gordon, position = position_brain(hemi ~ side), aes(fill = weights)) +
      scale_fill_gradientn(colors=rev(colfunc(16))) +
      facet_grid(rows=vars(dimensions),cols=vars(feature_type)) + 
      theme(strip.text = element_text(color = "white", face="bold")) + 
      theme(strip.background = element_rect(colour = "black", fill = "black")) +
      theme(legend.position = "none",axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggsave(plot=FIGURE6,filename=paste0(FIGURES,"/F6_",TODAY,".png"),device="png",width=14,height=8,units='in')

###### 
### Supplemental Figure 1: Prediction Accuracy
######

SF1<-fPERDICT[grep("FEAT-",fPERDICT$model),]
SF1$model<-factor(gsub("FEAT-","",SF1$model))
SF1<-as.data.frame(pivot_longer(SF1,2:4,names_to="dimension",values_to="prediction"))
SF1$dimension<-factor(gsub("_predict","",SF1$dimension)) 
df_scores<-df[,c(1,grep("_score",names(df)))] ; names(df_scores)[1]<-"sub"
df_scores<-as.data.frame(pivot_longer(df_scores,2:4,names_to="dimension",values_to="score"))
df_scores$dimension<-factor(gsub("_score","",df_scores$dimension))
SF1<-merge(SF1,df_scores,by=c("sub","dimension"))
SF1B<-SF1[which(SF1$model=="Multimod"),] 
SF1A<-SF1[-which(SF1$model=="Multimod"),] 
SF1A$model<-factor(SF1A$model,levels=c("CortThick","NetCohesion","OD","CortVolume","ReHo","ICVF","SubVolume","FCON","MD","SurfArea","ALFF","FA"))
FIGURE1SA<-ggplot() + 
  geom_abline(intercept=0,slope=1,,size=2,color="#000000",linetype="dashed") + 
  geom_point(data=SF1A, aes(x=score, y=prediction, color=dimension), size=0.75, alpha=0.3) +
  geom_smooth(data=SF1A, aes(x=score, y=prediction, group=dimension, color=dimension), method=lm, se=F, size=1.5, alpha=0.5) + 
  scale_color_manual(values=c("#0000FF","#ff0000","#00ce02")) +
  theme(strip.background = element_rect(colour = "black", fill = "black")) +
  theme(strip.text = element_text(color = "white", face="bold")) + 
  facet_wrap_equal(~model, scale="free",ncol=3) +
  xlab("Neural Properties") + ylab("Model Accuracy (Mean Squared Error)") + 
  theme(legend.position="none", panel.grid.minor=element_blank(), axis.title=element_blank()) 
FIGURE1SB<-ggplot() + 
  geom_abline(intercept=0,slope=1,,size=3,color="#000000",linetype="dashed") + 
  geom_point(data=SF1B, aes(x=score, y=prediction, color=dimension), size=1.5, alpha=0.25) +
  geom_smooth(data=SF1B, aes(x=score, y=prediction, group=dimension, color=dimension), fullrange=TRUE, method=lm, se=F, size=3, alpha=0.5) + 
  scale_color_manual(values=c("#0000FF","#ff0000","#00ce02")) +
  theme(strip.background = element_rect(colour = "black", fill = "black")) +
  theme(strip.text = element_text(color = "white", face="bold")) + 
  xlab("Psychopathology Symptoms") + ylab("Psychopathology Predictions") + 
  xlim(0,5.1) + ylim(0,5.1) + theme(legend.position="none") 
ggsave(plot=FIGURE1SB,filename=paste0(FIGURES,"/FS1_",TODAY,".png"),device="png",width=10.5,height=10.5,units='in')

###### 
### Supplemental Figure 2: Neurodevelopmental Trajectories
######

SF2<-df[,c(20,23:2341)]
for (INDEX in 2:2320){
  SF2[,INDEX]<-scale(SF2[,INDEX])
}
SF2<-pivot_longer(SF2,2:2320,names_to="features",values_to="neuro_score")
SF2<-merge(SF2,featmap,by="features") 
SF2$feature_type<-factor(SF2$feature_type)
SF2<-SF2[order(SF2$feature_type),]

FIGURE2S<-ggplot() +
  geom_smooth(data=SF2,aes(x=pfact_score,y=neuro_score,group=features,color=feature_type),method="loess",se=F,size=0.08) +
  scale_color_manual(values=unique(SF2$feature_color)) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour = "black", fill = "black")) +
  theme(strip.text = element_text(color = "white", face="bold")) + 
  xlab("General Psychopathology Factor (p-factor)") + ylab("Neuroimaging Metrics (Z-scored)") 
ggsave(plot=FIGURE2S,filename=paste0(FIGURES,"/FS2_",TODAY,".png"),device="png",width=10.5,height=10.5,units='in')

#####
### Supplemental Figure 3: Correlation Between Psychopathology Domains
#####

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
    } ; colnames(p.mat) <- rownames(p.mat) <- colnames(mat) ; p.mat
}

MAT<-df[,c("age","inter_score","exter_score","pfact_score")]
names(MAT)<-c("Age","Internalizing","Externalizing","p-factor")
r.mat<-as.matrix(cor(MAT,use="complete.obs")) ; p.mat <- cor.mtest(MAT)
col<-colorRampPalette(c("#1c4670","#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444","#8f1010"))
png(paste0(FIGURES,"/FS3_",TODAY,".png"),width=1008,height=1008,pointsize=24) 
corrplot(r.mat, method="color",type="lower",col=col(5000),addCoef.col="black",tl.col="black",tl.srt=45,p.mat=p.mat,sig.level=0.05,insig="blank",diag=T)
dev.off() 

#####
### Statistics for Tables
#####

#Model Accuracy
RESAMPLE<-read.csv(paste0(ROOT,"/Aggregate/RESAMPLE.csv"))
SUMMARIZE<-as.data.frame(RESAMPLE %>%
  group_by(feature_type, dimension) %>%
  summarize(
      mse_avg=mean(mse),low_mse=t.test(mse,conf.level = 0.99)$conf.int[1], high_mse=t.test(mse,conf.level = 0.95)$conf.int[2],
      mae_avg=mean(mae), low_mae=t.test(mae,conf.level = 0.99)$conf.int[1], high_mae=t.test(mae,conf.level = 0.95)$conf.int[2],
      corr_avg=mean(cod), low_cod=t.test(cod,conf.level = 0.99)$conf.int[1], high_cod=t.test(cod,conf.level = 0.95)$conf.int[2],
      .groups='drop'
    ) %>%
  mutate(across(where(is.numeric), round, digits = 3))
)

T1<-SUMMARIZE[which(SUMMARIZE$dimension=="Internalized Symptoms"),]
T1<-T1[grep("NET",T1$feature_type),] ; T1$dimension<-NULL
T1[order(T1$mse_avg),c(1,8:10)]

T2<-SUMMARIZE[which(SUMMARIZE$dimension=="Externalized Symptoms"),]
T2<-T2[grep("NET",T2$feature_type),] ; T2$dimension<-NULL
T2[order(T2$mse_avg),c(1,8:10)]

T3<-SUMMARIZE[which(SUMMARIZE$dimension=="General Psychopathology Factor"),]
T3<-T3[grep("FEAT",T3$feature_type),] ; T3$dimension<-NULL
T3[order(T3$mse_avg),c(1,8:10)]

RESAMPLE$dimension<-factor(RESAMPLE$dimension)
for (TYPE in unique(RESAMPLE$feature_type)){
  print(TYPE)
  SUBSET<-RESAMPLE[which(RESAMPLE$feature_type == TYPE),]
  print(eta_squared(lm(cod ~ dimension, data=SUBSET), partial = FALSE))
  #print(model_parameters(aov(cod ~ dimension, data=SUBSET)))
}
TEXT<-RESAMPLE[which(RESAMPLE$dimension=="General Psychopathology Factor"),]
TEXT<-TEXT[grep("NET",TEXT$feature_type),] ; TEXT<-TEXT[!grepl("NET-Global",TEXT$feature_type),]
model_parameters(aov(mse ~ feature_type, data=TEXT))
eta_squared(lm(mse ~ feature_type, data=TEXT))

#Gini Impurities
fEXPLAIN<-merge(fEXPLAIN,featmap[,c("features","feature_type")],by="features")
Rfeatmap<-featmap[!is.na(featmap$network_label),]
Rfeatmap$network_label<-gsub("NET-","",Rfeatmap$network_label)
nEXPLAIN<-merge(nEXPLAIN,Rfeatmap,by="features") 
fEXPLAIN$importance<-fEXPLAIN$importance*100
nEXPLAIN$importance<-nEXPLAIN$importance*100
SUMMARIZE<-as.data.frame(nEXPLAIN %>%
  group_by(network_label, dimensions) %>%
  summarize(
      #avg=mean(importance), median=median(importance), sd=sd(importance), sum=sum(importance),
      avg=mean(importance), low=t.test(importance,conf.level = 0.99)$conf.int[1], high=t.test(importance,conf.level = 0.99)$conf.int[2],
      .groups='drop'
    ) 
) 
SUMMARIZE$median<-round(SUMMARIZE$median, digits=3)
SUMMARIZE$avg<-round(SUMMARIZE$avg, digits=3)
SUMMARIZE$sd<-round(SUMMARIZE$sd, digits=2)
SUM<-SUMMARIZE[which(SUMMARIZE$dimension=="pfact_score"),]
SUM[order(SUM$median),]

for (TYPE in levels(factor(fEXPLAIN$feature_type))){
  SUBSET<-fEXPLAIN[which(fEXPLAIN$feature_type==TYPE),]
  print(TYPE)
  print(eta_squared(lm(importance ~ dimensions, data=SUBSET)))
}

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######