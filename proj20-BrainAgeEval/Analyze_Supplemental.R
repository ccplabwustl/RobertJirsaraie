#!/usr/bin/env Rscript
######################

library(gamm4) ; library(lme4) ; library(nlme) ; library(parameters) ; library(performance)
library(ggplot2) ; library(ggforce) ; library(corrplot) ; library(ggeffects) 
library(tidyr) ; library(plyr) ; library(dplyr) ; library(purrr) 
library(knitr) ; library(psych) ; library(broom.mixed) 
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
### Read the Datasets
######

NEURO_PATH<-"/Users/Jirsaraie/Box/Research/proj20-BrainAgeEval/audits"
LONG<-read.csv(paste0(NEURO_PATH,"/NEURO_Traj.csv"))
FINAL<-read.csv(paste0(NEURO_PATH,"/NEURO_Grp.csv"))
COMBINE<-read.csv(paste0(NEURO_PATH,"/COMBINE.csv"))
PRISMA<-read.csv(paste0(NEURO_PATH,"/PRISMA.csv"))
TRIO<-read.csv(paste0(NEURO_PATH,"/TRIO.csv"))
HCP<-read.csv(paste0(NEURO_PATH,"/HCP.csv"))

###Scatterplots of Brain Age Models Fit with Age By Dataset

HCP<-HCP[,-c(69:74)]
for (COL in 4:9){
	ERROR<-paste0(names(TRIO)[COL],"_AgeAdjustedError")
	TRIO[,ERROR]<-residuals(lmer(TRIO[,COL]~TRIO[,"age"]+ (1 | TRIO[,"sub"])))
	PRISMA[,ERROR]<-residuals(lmer(PRISMA[,COL]~PRISMA[,"age"]+ (1 | PRISMA[,"sub"])))
	HCP[complete.cases(HCP[,COL]),ERROR]<-lm(HCP[,COL]~HCP[,"age"])$resid
	PERDICTION<-paste0(names(TRIO)[COL],"_AgeAdjustedPred")
	TRIO[,PERDICTION]<-residuals(lmer(TRIO[,COL]~TRIO[,"age"]+ (1 | TRIO[,"sub"])))+TRIO[,"age"]
	PRISMA[,PERDICTION]<-residuals(lmer(PRISMA[,COL]~PRISMA[,"age"]+ (1 | PRISMA[,"sub"])))+PRISMA[,"age"]
	HCP[complete.cases(HCP[,COL]),PERDICTION]<-lm(HCP[,COL]~HCP[,"age"])$resid+HCP[complete.cases(HCP[,COL]),"age"]
	NORM<-paste0(names(TRIO)[COL],"_AgeAdjustedNorm")
	TRIO[,NORM]<-normalize(residuals(lmer(TRIO[,COL]~TRIO[,"age"]+ (1 | TRIO[,"sub"]))))
	PRISMA[,NORM]<-normalize(residuals(lmer(PRISMA[,COL]~PRISMA[,"age"]+ (1 | PRISMA[,"sub"]))))
	HCP[complete.cases(HCP[,COL]),NORM]<-normalize(lm(HCP[,COL]~HCP[,"age"])$resid)
}
COMBINE<-rbind.fill(TRIO,PRISMA,HCP) 

###Scatterplots of Brain Age Models Fit with Age By Dataset
LONG <-pivot_longer(COMBINE,cols=grep("AdjustedPred",names(COMBINE)),names_to="MODELS",values_to="BRAINAGE") 
LONG$MODELS<-factor(gsub("_AgeAdjustedPred","",LONG$MODELS),levels=c("DBN","rDBN","tDBN","GTB","rGTB","tGTB"))
ggplot(LONG, aes(x=age,y=BRAINAGE,group=DATA,color=DATA,shape=DATA)) + 
	geom_abline(intercept=0,slope=1,size=2,alpha=1,linetype="dashed") +
	geom_point(size=1.45,alpha=0.575) +
	stat_smooth(method = "lm",se=F,size=2.5,alpha=0.5,fullrange=T) +
	scale_shape_manual(values=c(16,17,15))+
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS,scales = "free")

###Benchmark Developmental Brain Age Metrics
for (COL in grep("AdjustedPred",names(COMBINE))){
	DATA<-TRIO ; TYPE<-"LONG" ; print(paste0(names(DATA)[COL]))
	if (TYPE == "CROSS"){
		COEF<-round(lm(DATA[,COL]~DATA[,"age"])$coefficients,digits=2)
	} else if (TYPE == "LONG"){
		COEF<-suppressMessages(round(fixef(lmer(DATA[,COL] ~ DATA[,"age"] + (1 | DATA[,"sub"]))),digits=2))
	}
	MAE<-round(mean(abs((DATA[,COL]-DATA[,"age"])),na.rm=T),digits=2)
	CORR<-round(cor(DATA[,COL],DATA[,"age"],use = "complete.obs"),digits=2)
	print(paste0("MAEx",MAE,"_SLOPEx",COEF[2],"_Interceptx",COEF[1],"_CORRx",CORR))
}

###Add Reliability Metrics Based on Corrected Brain Ages
for (ROW in 1:dim(COMBINE)[1]){
	COMBINE[ROW,"DBN_SD_AgeAdj_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_AgeAdjustedError",names(COMBINE))[1:3]]]),na.rm=T)
	COMBINE[ROW,"DBN_SD_AgeAdj_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_AgeAdjustedNorm",names(COMBINE))][1:3]]),na.rm=T)
	COMBINE[ROW,"GTB_SD_AgeAdj_RAW"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_AgeAdjustedError",names(COMBINE))[4:6]]]),na.rm=T)
	COMBINE[ROW,"GTB_SD_AgeAdj_NORM"]<-sd((COMBINE[ROW,names(COMBINE)[grep("_AgeAdjustedNorm",names(COMBINE))][4:6]]),na.rm=T)
}
COMBINE$DATA<-factor(COMBINE$DATA,levels=c("PDS_TRIO","PDS_PRISMA","HCP_TEST"))
PRISMA<-COMBINE[which(COMBINE$DATA=="PDS_PRISMA"),] 
TRIO<-COMBINE[which(COMBINE$DATA=="PDS_TRIO"),] 
HCP<-COMBINE[which(COMBINE$DATA=="HCP_TEST"),]

#Histogram of Reliability by Models
HIST_RAW <-pivot_longer(COMBINE,cols=c(87,89),names_to="MODELS",values_to="DEVS") 
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_SD_AgeAdj_RAW), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_SD_AgeAdj_RAW), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_RAW, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")
HIST_NORM <-pivot_longer(COMBINE,cols=c(88,90),names_to="MODELS",values_to="DEVS")
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_SD_AgeAdj_NORM), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_SD_AgeAdj_NORM), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_NORM, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")

###Get Model Parameters Regarding Individual Differences in Relability
model_parameters(lmer(DBN_SD_AgeAdj_NORM ~ age + sex + EulerNumber + (1 | sub),data=PRISMA))
model_parameters(lmer(DBN_SD_AgeAdj_NORM ~ age + sex + EulerNumber +  (1 | sub),data=TRIO))
model_parameters(lm(DBN_SD_AgeAdj_NORM ~ age + sex + EulerNumber + numNavs_T1w, data=HCP)) 
LONG <-pivot_longer(COMBINE,cols=grep("_SD_AgeAdj_NORM",names(COMBINE)),names_to="MODELS",values_to="DEVS")
ggplot(LONG, aes(x=age,y=DEVS,group=DATA,color=DATA,shape=DATA)) + 
	geom_point(size=1.5,alpha=0.7) + 	scale_shape_manual(values=c(16,17,15)) +
	stat_smooth(method = "loess",se=F,size=3,alpha=1,fullrange=T) +
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS, scale="free")

###Get Standardized Betas
LmerBetas(lmer(DBN_SD_AgeAdj_NORM ~ age + EulerNumber + sex + (1|sub), data = TRIO))
LmerBetas(lmer(DBN_SD_AgeAdj_NORM ~ age + EulerNumber + sex + (1|sub), data = PRISMA))
lm.beta(lm(DBN_SD_AgeAdj_NORM ~ age + EulerNumber + sex + numNavs_T1w, data = HCP))

LmerBetas(lmer(GTB_SD_AgeAdj_NORM ~ age + EulerNumber + sex + (1|sub), data = TRIO))
LmerBetas(lmer(GTB_SD_AgeAdj_NORM ~ age + EulerNumber + sex + (1|sub), data = PRISMA))
lm.beta(lm(GTB_SD_AgeAdj_NORM ~ age + EulerNumber + sex + numNavs_T1w, data = HCP))

LmerBetas(lmer(DBN_SD_AgeAdj_NORM ~ DBN_SD_NORM + (1|sub), data = TRIO))
LmerBetas(lmer(DBN_SD_AgeAdj_NORM ~ DBN_SD_NORM + (1|sub), data = PRISMA))
lm.beta(lm(DBN_SD_AgeAdj_NORM ~ DBN_SD_NORM, data = HCP))

###Non-Linearity Tests
exactRLRT(gamm(DBN_SD_NORM~s(age), data=HCP)$lme)

###Add Reliability Metrics Based on Skewness
for (ROW in 1:dim(COMBINE)[1]){
	COMBINE[ROW,"DBN_SKEW_RAW"]<-skewness(t(COMBINE[ROW,names(COMBINE)[grep("_ERROR",names(COMBINE))[1:3]]]))[1,2]
	COMBINE[ROW,"DBN_SKEW_NORM"]<-skewness(t(COMBINE[ROW,names(COMBINE)[grep("_ERROR_NORM",names(COMBINE))][1:3]]))[1,2]
	COMBINE[ROW,"GTB_SKEW_RAW"]<-skewness(t(COMBINE[ROW,names(COMBINE)[grep("_ERROR",names(COMBINE))[4:6]]]))[1,2]
	COMBINE[ROW,"GTB_SKEW_NORM"]<-skewness(t(COMBINE[ROW,names(COMBINE)[grep("_ERROR_NORM",names(COMBINE))][4:6]]))[1,2]
}
COMBINE$DATA<-factor(COMBINE$DATA,levels=c("PDS_TRIO","PDS_PRISMA","HCP_TEST"))
PRISMA<-COMBINE[which(COMBINE$DATA=="PDS_PRISMA"),] 
TRIO<-COMBINE[which(COMBINE$DATA=="PDS_TRIO"),] 
HCP<-COMBINE[which(COMBINE$DATA=="HCP_TEST"),]

###Make Histograms of Skewness
HIST_RAW <-pivot_longer(COMBINE,cols=c(91,93),names_to="MODELS",values_to="DEVS") 
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_SKEW_RAW,na.rm=T), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_SKEW_RAW,na.rm=T), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_RAW, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")
HIST_NORM <-pivot_longer(COMBINE,cols=c(92,94),names_to="MODELS",values_to="DEVS")
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_SKEW_NORM,na.rm=T), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_SKEW_NORM,na.rm=T), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_NORM, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")
LONG <-pivot_longer(COMBINE,cols=grep("_SKEW_NORM",names(COMBINE)),names_to="MODELS",values_to="DEVS")
ggplot(LONG, aes(x=age,y=DEVS,group=DATA,color=DATA,shape=DATA)) + 
	geom_point(size=1.5,alpha=0.7) + 	scale_shape_manual(values=c(16,17,15)) +
	stat_smooth(method = "loess",se=F,size=3,alpha=1,fullrange=T) +
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS, scale="free")

###Find Correlation Between Raw and Corrected Brain Ages
LONG <-pivot_longer(COMBINE,cols=4:9,names_to="MODELS",values_to="BRAINAGE") 
LONG$MODELS<-factor(LONG$MODELS,levels=c("DBN","rDBN","tDBN","GTB","rGTB","tGTB"))
for(MOD in levels(LONG$MODELS)){
	ADJUSTEDCOL<-which(names(LONG)==paste0(MOD,"_AgeAdjustedPred"))
	VALUES<-LONG[which(LONG$MODELS==MOD),ADJUSTEDCOL]
	LONG[which(LONG$MODELS==MOD),"CORRECTAGE"]<-VALUES
}
LONG$DATA<-factor(LONG$DATA,levels=c("PDS_TRIO","PDS_PRISMA","HCP_TEST"))
lPRISMA<-LONG[which(LONG$DATA=="PDS_PRISMA"),] 
lTRIO<-LONG[which(LONG$DATA=="PDS_TRIO"),] 
lHCP<-LONG[which(LONG$DATA=="HCP_TEST"),]
for(MOD in levels(lHCP$MODELS)){
	TEMP<-lHCP[which(lHCP$MODELS == MOD),]
	COR<-cor(TEMP$BRAINAGE,TEMP$CORRECTAGE,use = "complete.obs")
	print(paste0(MOD,"x",COR))
}
ggplot(LONG, aes(x=CORRECTAGE,y=BRAINAGE,group=DATA,color=DATA,shape=DATA)) + 
	geom_abline(intercept=0,slope=1,size=2,alpha=1,linetype="dashed") +
	geom_point(size=1.45,alpha=0.575) +
	stat_smooth(method = "lm",se=F,size=2.5,alpha=0.5,fullrange=T) +
	scale_shape_manual(values=c(16,17,15))+
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS,scales = "free")

###Recreate Figure 2 with Ranked Brain Ages Instead of BAGs
LONG <-pivot_longer(COMBINE,cols=4:9,names_to="MODELS",values_to="BRAINAGE") #4:9 - 32:37 - 38:43
ggplot(LONG, aes(x=age, y=BRAINAGE, color=MODELS)) + geom_point(size=1.75) + theme_classic() +
	scale_color_manual(values=c("#0022ff","#00d9ff","#040080","#7c3f00","#db9e65","#472309")) +
	facet_wrap(~LONG$DATA,scale="free")
for (COL in 4:9){
	TRIO[,paste0(names(TRIO)[COL],"_NORM")]<-normalize(COMBINE[,COL])
	COMBINE[,paste0(names(COMBINE)[COL],"_NORM")]<-normalize(COMBINE[,COL])
	COMBINE[,paste0(names(COMBINE)[COL],"_NORM")]<-normalize(COMBINE[,COL])
}
LONG <-pivot_longer(COMBINE,cols=(ncol(COMBINE)-5):ncol(COMBINE),names_to="MODELS",values_to="BRAINAGE") #4:9 - 32:37 - 38:43
ggplot(LONG, aes(x=age, y=BRAINAGE, color=MODELS)) + geom_point(size=1.75) + theme_classic() +
	scale_color_manual(values=c("#0022ff","#00d9ff","#040080","#7c3f00","#db9e65","#472309")) +
	facet_wrap(~LONG$DATA,scale="free")

###Add Reliability Metrics Based on Coefficient of Variation
for (ROW in 1:dim(COMBINE)[1]){
	COMBINE[ROW,"DBN_CoefVar_RAW"]<-sd(COMBINE[ROW,grep("_ERROR",names(COMBINE))[1:3]])/mean(t(COMBINE[ROW,grep("_ERROR",names(COMBINE))[1:3]]))
	COMBINE[ROW,"GTB_CoefVar_RAW"]<-sd(COMBINE[ROW,grep("_ERROR",names(COMBINE))[4:6]])/mean(t(COMBINE[ROW,grep("_ERROR",names(COMBINE))[4:6]]))
	COMBINE[ROW,"DBN_CoefVar_NORM"]<-sd(COMBINE[ROW,grep("_ERROR",names(COMBINE))[7:9]])/mean(t(COMBINE[ROW,grep("_ERROR",names(COMBINE))[7:9]]))
	COMBINE[ROW,"GTB_CoefVar_NORM"]<-sd(COMBINE[ROW,grep("_ERROR",names(COMBINE))[10:12]])/mean(t(COMBINE[ROW,grep("_ERROR",names(COMBINE))[10:12]]))
}
COMBINE$DATA<-factor(COMBINE$DATA,levels=c("PDS_TRIO","PDS_PRISMA","HCP_TEST"))
PRISMA<-COMBINE[which(COMBINE$DATA=="PDS_PRISMA"),] 
TRIO<-COMBINE[which(COMBINE$DATA=="PDS_TRIO"),] 
HCP<-COMBINE[which(COMBINE$DATA=="HCP_TEST"),]
HIST_RAW <-pivot_longer(COMBINE,cols=grep("CoefVar_RAW",names(COMBINE)),names_to="MODELS",values_to="DEVS") 
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_CoefVar_RAW,na.rm=T), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_CoefVar_RAW,na.rm=T), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_RAW, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")
HIST_NORM <-pivot_longer(COMBINE,cols=grep("CoefVar_NORM",names(COMBINE)),names_to="MODELS",values_to="DEVS")
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_CoefVar_NORM,na.rm=T), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_CoefVar_NORM,na.rm=T), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_NORM, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")
LONG <-pivot_longer(COMBINE,cols=grep("CoefVar_NORM",names(COMBINE)),names_to="MODELS",values_to="DEVS")
ggplot(LONG, aes(x=age,y=DEVS,group=DATA,color=DATA,shape=DATA)) + 
	geom_point(size=1.5,alpha=0.7) + 	scale_shape_manual(values=c(16,17,15)) +
	stat_smooth(method = "loess",se=F,size=3,alpha=1,fullrange=T) +
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS, scale="free")
###Visualize the Age-Dependency Artifact for all Models
model_parameters(gamm4(DBN_CoefVar_NORM ~ s(age) + sex + EulerNumber, random=~(1 | sub), data=PRISMA))
model_parameters(gamm4(DBN_CoefVar_NORM ~ s(age) + sex + EulerNumber, random=~(1 | sub), data=TRIO))
model_parameters(gam(DBN_CoefVar_NORM ~ s(age) + sex + EulerNumber, data=HCP))
model_parameters(gamm4(GTB_CoefVar_NORM ~ s(age) + sex + EulerNumber, random=~(1 | sub), data=PRISMA))
model_parameters(gamm4(GTB_CoefVar_NORM ~ s(age) + sex + EulerNumber, random=~(1 | sub), data=TRIO))
model_parameters(gam(GTB_CoefVar_NORM ~ s(age) + sex + EulerNumber, data=HCP))

###Reliability Differences Between Models
DIFF<-pivot_longer(HCP,cols=grep("CoefVar_NORM",names(HCP)),names_to="MODELS",values_to="DEVS")
t.test(DIFF$DEVS~DIFF$MODELS)

###Create Correlation Matrix Between All Brain Age Models
MAT<-HCP[,c(2:8)] ; MAT<-TRIO[,c(3:9)] ; MAT<-PRISMA[,c(3:9)] 
names(MAT)[1]<-"Age" ; r.mat<-as.matrix(cor(MAT,use="complete.obs")) ; p.mat <- cor.mtest(MAT) 
col<-colorRampPalette(c("#1C4670","#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444","#8F1010"))
corrplot(r.mat,method="color",type="lower",col=col(10000),addCoef.col="#FFFFFF",tl.col="black", tl.srt=45,sig.level=0.01,insig="blank",diag=T)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######