#!/usr/bin/env Rscript
######################

TODAY<-format(Sys.time(), "%Y%m%d")
library(stringr) ; library(dplyr) ; library(tidyverse) ; library(corrplot) ; library(lattice) ; library(lme4)
COLLAB_DIR<-"/Users/Jirsaraie/Desktop/Research/collaborations/20A-KLuking/Data"
DATA<-read.csv(list.files(COLLAB_DIR,full.names=T,pattern=".csv")[1])
SCAN1<-read.csv(list.files(COLLAB_DIR,full.names=T,pattern=".csv")[2]) 
SCAN2<-read.csv(list.files(COLLAB_DIR,full.names=T,pattern=".csv")[3]) 
SCAN3<-read.csv(list.files(COLLAB_DIR,full.names=T,pattern=".csv")[4])
SCAN1$ses<-1 ; names(SCAN1)<-gsub("PDS_SCAN1_Penn_","",names(SCAN1))
SCAN2$ses<-2 ; names(SCAN2)<-gsub("PDS_SCAN2_Penn_","",names(SCAN2))
SCAN3$ses<-3 ; names(SCAN3)<-gsub("PDS_SCAN3_Penn_","",names(SCAN3))
LONG<-rbind(SCAN1,SCAN2,SCAN3) ; names(LONG)[1]<-c("sub") ; LONG<-LONG[,c(1,20,2:19)]
LONG$sub<-as.numeric(gsub("B","",gsub("L","",LONG$sub)))
MASTER<-read.csv("/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/analysis/n834_DataFreeze_20201231.csv")
MASTER<-MASTER[,c("sub","ses","age")] ; LONG<-merge(MASTER,LONG,by=c("sub","ses"),all=F)
DATA[,1:ncol(DATA)]<-lapply(DATA[,1:ncol(DATA)],as.numeric)
DATA$sex<-as.factor(DATA$sex)

######
### Create Figure 1 - Timevarying Cortical Thickness Trajectories
######
#solid dashed dotted dotdash dotted
#000000         9,10,14,17,18
#006393 blue    1,5,13,16
#cb5700 orange  8
#4da1c5 lblue   7,11,12 
#048c61 green   6
#d5c637 yellow  2,3,15
#b36a8a pink    4

FIGURE1 <- LONG %>% pivot_longer(cols = 4:21, names_to = "Component", values_to = "Thickness")
COLORS <- c("#006393","#d5c637","#d5c637","#b36a8a","#006393","#048c61","#4da1c5","#cb5700","#000000","#000000","#4da1c5","#4da1c5","#006393","#000000","#d5c637","#006393","#000000","#000000")
SHAPE <- c("twodash","longdash","dashed","solid","longdash","solid","longdash","solid","solid","twodash","dashed","solid","dashed","longdash","solid","solid","dashed","dotted")

#linetype=Component
ggplot(FIGURE1, aes(x=age,y=Thickness,group=Component,color=Component)) + 
	geom_smooth(method = lm,se=F,size=1.25,alpha=5) +
	scale_linetype_manual(values=SHAPE) + 
	scale_color_manual(values=COLORS) + 
	theme_classic()

######
### Create Figure 2 - ScatterPlots of NeuroBehav Relationships
######

ANX<-as.data.frame(matrix(NA,ncol=6,nrow=14)) ; INDEX=0
for (COMPNUM in c(2,3,4,6,15)){
	INDEX=((INDEX+1))
	COLNUM<-grep(paste0("Cmp_",COMPNUM,"_int"),names(DATA))
	MODEL<-lm(DATA[,COLNUM]~S1AGEMO+sex+T1Income_to_Need+AVG_MDD_C_PS+AVG_MDD_C_SA+AVG_MDD_C_AD+EXTLdim_PS+EXTLdim_SA+EXTLdim_AD+INTLdim_PS+INTLdim_SA+INTLdim_AD,data=DATA)
	OUTPUT<-ggpredict(MODEL, terms = "INTLdim_PS")
	ANX[,INDEX]<-OUTPUT$predicted ; names(ANX)[INDEX]<-paste0(names(DATA)[COLNUM],"_ANX")
}
ANX[,INDEX+1]<-OUTPUT$x ; names(ANX)[INDEX+1]<-"X" ; ANX$SCATTER<-1
ANX<-ANX %>% pivot_longer(cols = 1:5, names_to = "relationship", values_to = "value")

POV<-as.data.frame(matrix(NA,ncol=7,nrow=9)) ; INDEX=0
for (COMPNUM in c(4,6,7,8,11,12)){
	INDEX=((INDEX+1))
	COLNUM<-grep(paste0("Cmp_",COMPNUM,"_int"),names(DATA))
	MODEL<-lm(DATA[,COLNUM]~S1AGEMO+sex+T1Income_to_Need+AVG_MDD_C_PS+AVG_MDD_C_SA+AVG_MDD_C_AD+EXTLdim_PS+EXTLdim_SA+EXTLdim_AD+INTLdim_PS+INTLdim_SA+INTLdim_AD,data=DATA)
	OUTPUT<-ggpredict(MODEL, terms = "T1Income_to_Need")
	POV[,INDEX]<-OUTPUT$predicted ; names(POV)[INDEX]<-paste0(names(DATA)[COLNUM],"_POV")
}
POV[,INDEX+1]<-OUTPUT$x ; names(POV)[INDEX+1]<-"X" ; POV$SCATTER<-2
POV<-POV %>% pivot_longer(cols = 1:6, names_to = "relationship", values_to = "value")


EXT<-as.data.frame(matrix(NA,ncol=7,nrow=8)) ; INDEX=0
for (COMPNUM in c(1,5,6,8,13,16)){
	INDEX=((INDEX+1))
	COLNUM<-grep(paste0("Cmp_",COMPNUM,"_slp"),names(DATA))
	MODEL<-lm(DATA[,COLNUM]~S1AGEMO+S3AGEMO+sex+T1Income_to_Need+AVG_MDD_C_PS+AVG_MDD_C_SA+AVG_MDD_C_AD+EXTLdim_PS+EXTLdim_SA+EXTLdim_AD+INTLdim_PS+INTLdim_SA+INTLdim_AD,data=DATA)
	OUTPUT<-ggpredict(MODEL, terms = "EXTLdim_PS")
	EXT[,INDEX]<-OUTPUT$predicted ; names(EXT)[INDEX]<-paste0(names(DATA)[COLNUM],"_EXT")
}
EXT[,INDEX+1]<-OUTPUT$x ; names(EXT)[INDEX+1]<-"X" ; EXT$SCATTER<-3
EXT<-EXT %>% pivot_longer(cols = 1:6, names_to = "relationship", values_to = "value")
FIGURE2<-rbind(ANX,POV,EXT)


COLORS <- c("#d5c637","#d5c637","#b36a8a","#048c61","#d5c637","#b36a8a","#048c61","#4da1c5","#cb5700","#4da1c5","#4da1c5","#006393","#006393","#048c61","#cb5700","#006393","#d5c637")
SHAPE <- c("longdash","dashed","solid","solid","solid","solid","longdash","solid","solid","twodash","dashed","solid","dashed","longdash","solid","solid","dashed")

#linetype=relationship
ggplot(FIGURE2, aes(x=X,y=value,group=relationship,color=relationship)) + 
	geom_smooth(method = lm,se=F,size=1.25,alpha=5) +
	scale_linetype_manual(values=SHAPE) + 
	scale_color_manual(values=COLORS) + 
	facet_wrap(~SCATTER, scales = "free") + 
	theme_classic() 

######
### Create Figure 3 - Correlation Matrix of Predictors
######

names(DATA)[which(names(DATA)=="Subid")]<-"sub"
TEMP1<-LONG[which(LONG$ses == 1),c("sub","age")] ; names(TEMP1)[2]<-"age1"
TEMP2<-LONG[which(LONG$ses == 2),c("sub","age")] ; names(TEMP2)[2]<-"age2"
TEMP3<-LONG[which(LONG$ses == 3),c("sub","age")] ; names(TEMP3)[2]<-"age3"
DATA<-merge(DATA,TEMP1,by=c("sub")) ; DATA<-merge(DATA,TEMP2,by=c("sub")) ; DATA<-merge(DATA,TEMP3,by=c("sub"))
MAT<-DATA[,c("age1","age2","age3","T1Income_to_Need","AVG_MDD_C_PS","AVG_MDD_C_SA","AVG_MDD_C_AD","EXTLdim_PS","EXTLdim_SA","EXTLdim_AD","INTLdim_PS","INTLdim_SA","INTLdim_AD")]
names(MAT)<-c("ScanAge1","ScanAge2","ScanAge3","IncomeToNeeds","MDD_PS","MDD_SA","MDD_AD","EXTL_PS","EXTL_SA","EXTL_AD","INTL_PS","INTL_SA","INTL_AD")

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

r.mat<-as.matrix(cor(MAT,use="complete.obs")) ; p.mat <- cor.mtest(MAT)
col <- colorRampPalette(c("#8f1010","#BB4444","#EE9988","#FFFFFF","#77AADD","#4477AA","#1c4670"))
col <- colorRampPalette(c("#1c4670","#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444","#8f1010"))
corrplot(r.mat, method="color",type="lower",col=col(5000),addCoef.col="black",tl.col="black", tl.srt=45, p.mat = p.mat, sig.level = 0.05, insig = "blank",diag=T)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######