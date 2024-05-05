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

######
### Analyze ANOVAs Using Groups of Brain Age from the HCP Dataset
######

###Alternative Approach with Residuals
LONG <- pivot_longer(HCP, cols = 3:8, names_to="models", values_to="brainages") 
LONG_NEST <- LONG %>% group_by(models) %>% nest() ; row.names(HCP)<-NULL
MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lm(brainages ~ age, data = .))) 
MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))
for (MOD in 1:dim(MODELS)[1]){
	RESID<-paste0(unique(LONG$models)[MOD],"_RESID") ; LABEL<-paste0(unique(LONG$models)[MOD],"_GROUP")
	DATASET<-as.data.frame(resid(MODELS$fit1[[MOD]])) ; names(DATASET)[1]<-RESID
	DATASET$COUNT<-as.numeric(row.names(DATASET))
	HCP$COUNT<-as.numeric(row.names(HCP)) 
	HCP<-merge(HCP,DATASET,by="COUNT",all=TRUE)
	HCP[which(HCP[,RESID] <= sort(HCP[,RESID])[round(nrow(DATASET)*0.30)]),LABEL]<-"Delayed"
	HCP[which(abs(HCP[,RESID]) <= sort(abs(HCP[,RESID]))[round(nrow(DATASET)*0.30)]),LABEL]<-"Normal"
	HCP[which(HCP[,RESID] >= sort(HCP[,RESID])[(nrow(DATASET)-round(nrow(DATASET)*0.30)+1)]),LABEL]<-"Accelerated"
	HCP[which(is.na(HCP[,LABEL])==TRUE),LABEL]<-"Exclude" ; HCP$COUNT<-NULL
	HCP[,LABEL] <- factor(HCP[,LABEL],levels=c("Delayed","Normal","Accelerated","Exclude")) 
}
LONG <- as.data.frame(pivot_longer(HCP, cols = grep("_GROUP",names(HCP)), names_to="MODELS", values_to="GROUPS"))
for (SUB in unique(LONG$sub)){
	NUMBER<-length(unique(LONG[which(LONG$sub==SUB),"GROUPS"]))
	if (NUMBER==1){
		HCP[which(HCP$sub==SUB),"GRP_CHG"]<-HCP[which(HCP$sub==SUB),"GTB_GROUP"]
	}
}

##### Visualize Developmental Trajectories that Parcelated Groups 
LONG <- HCP %>% pivot_longer(cols = 3:8, names_to="models", values_to="brainage")
LONG <- LONG %>% pivot_longer(cols = grep("_GROUP",names(LONG)), names_to="developmental", values_to="group")  
LONG$developmental<-gsub("_GROUP","",LONG$developmental) ; LONG<-LONG[which(LONG$models == LONG$developmental),] 
LONG$models<-factor(LONG$models, levels = c("DBN","rDBN","tDBN","GTB","rGTB","tGTB")) 
ggplot(LONG, aes(x=age,y=brainage)) + 
	geom_smooth(method = "lm",se=F,size=2,alpha=1,fullrange=F,color="#8B4513") +
	geom_point(size=2.5,alpha=0.75,aes(color=group)) +
	geom_abline(intercept=0,slope=1,size=2,alpha=1) +
	#coord_cartesian(xlim=c(5,25),ylim=c(5,25)) +
	scale_color_manual(values=c("#228B22", "#8B4513", "#e60000", "#b8b8b8")) + 
	theme_classic() + facet_wrap(~models,scale="free")

##### Visualize Developmental Trajectories that Parcelated Groups 
OUT<-as.data.frame(names(HCP)[3:8]) ; names(OUT)<-"models" ; LONG$models<-gsub("_GROUP","",LONG$models) 
LONG <-pivot_longer(HCP,cols=grep("_GROUP",names(HCP)),names_to="models",values_to="braingrps") 
LONG$braingrp<-LONG$braingrps ; LONG[which(LONG$braingrp=="Normal"),"braingrp"]<-NA 
LONG$braingrp<-factor(LONG$braingrp,levels=c("Delayed","Accelerated"))
LONG_NEST <- LONG %>% group_by(models) %>% nest
for (COL in c(2,19,21,20,18,17)){
	MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lm(HCP[,COL] ~ braingrp + EulerNumber + sex, data = .)))
	MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))  
	for (MOD in 1:dim(MODELS)[1]){
		print(paste0(MODELS$models[MOD],"x",names(HCP)[COL])) 
		print(model_parameters(MODELS$fit1[[MOD]]))
		OUT[which(OUT$models==gsub("_GROUP","",MODELS$models[MOD])),names(HCP)[COL]]<-round(model_parameters(MODELS$fit1[[MOD]])[2,6],digits=3) 
	}
}

######
### Analyze ANOVAs Using Groups of Brain Age from the PDS Dataset
######

LONG <- pivot_longer(PRISMA, cols = 4:9, names_to="models", values_to="brainages") 
LONG_NEST <- LONG %>% group_by(models) %>% nest() 
#MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lm(brainages ~ age, data = .))) 
MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lme(brainages ~ age, random = ~ 1 | sub, data = .)))
MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))
for (MOD in 1:dim(MODELS)[1]){
	RESID<-paste0(unique(LONG$models)[MOD],"_RESID") ; LABEL<-paste0(unique(LONG$models)[MOD],"_GROUP")
	DATASET<-as.data.frame(residuals(MODELS$fit1[[MOD]],level=0)) ; names(DATASET)[1]<-RESID
	DATASET$COUNT<-as.numeric(row.names(DATASET)) ; row.names(PRISMA)<-NULL
	PRISMA$COUNT<-as.numeric(row.names(PRISMA))
	PRISMA<-merge(PRISMA,DATASET,by="COUNT",all=TRUE)
	PRISMA[which(PRISMA[,RESID] <= sort(PRISMA[,RESID])[round(nrow(DATASET)*0.30)]),LABEL]<-"Delayed"
	PRISMA[which(abs(PRISMA[,RESID]) <= sort(abs(PRISMA[,RESID]))[round(nrow(DATASET)*0.30)]),LABEL]<-"Normal"
	PRISMA[which(PRISMA[,RESID] >= sort(PRISMA[,RESID])[(nrow(DATASET)-round(nrow(DATASET)*0.30)+1)]),LABEL]<-"Accelerated"
	PRISMA[which(is.na(PRISMA[,LABEL])==TRUE),LABEL]<-"Exclude" ; PRISMA$COUNT<-NULL
	PRISMA[,LABEL]<-factor(PRISMA[,LABEL],levels=c("Delayed","Normal","Accelerated","Exclude")) 
}

### Cluster Subjects into a Time-Invarying Developmental Group
for (SUB in unique(PRISMA$sub)){
	for (COL in grep("_GROUP",names(PRISMA))){
		LABEL<-gsub("_GROUP","_CLUSTER",names(PRISMA)[COL])
		if (length(unique(PRISMA[which(PRISMA$sub==SUB),COL]))==1){
			PRISMA[which(PRISMA$sub==SUB),LABEL]<-as.character(PRISMA[which(PRISMA$sub==SUB)[1],COL])
		} else {
			PRISMA[which(PRISMA$sub==SUB),LABEL]<-NA
		}
	}
}
PRISMA[,grep("CLUSTER",names(PRISMA))]<-lapply(PRISMA[,grep("CLUSTER",names(PRISMA))],factor)
LONG <- as.data.frame(pivot_longer(PRISMA, cols = grep("_CLUSTER",names(PRISMA))[1:6], names_to="MODELS", values_to="GROUPS"))
for (SUB in unique(PRISMA$sub)){
	NUMBER<-length(unique(LONG[which(LONG$sub==SUB),"GROUPS"]))
	if (NUMBER==1){
		PRISMA[which(PRISMA$sub==SUB),"GRP_CHG"]<-PRISMA[which(PRISMA$sub==SUB),"GTB_GROUP"]
	}
}

##### Visualize Developmental Trajectories that Parcelated Groups 
LONG <- PRISMA %>% pivot_longer(cols = 4:11, names_to="models", values_to="brainage")
LONG <- LONG %>% pivot_longer(cols = grep("_GROUP",names(LONG)), names_to="developmental", values_to="group")  
LONG$developmental<-gsub("_GROUP","",LONG$developmental) ; LONG<-LONG[which(LONG$models==LONG$developmental),] 
LONG$models<-factor(LONG$models, levels = c("DBN","rDBN","tDBN","GTB","rGTB","tGTB","hGTBgam","hGTBlmer")) 
ggplot(LONG, aes(x=age,y=brainage)) + 
	geom_smooth(method = "lm",se=F,size=2,alpha=1,fullrange=F,color="#8B4513") +
	geom_point(size=2.5,alpha=0.75,aes(color=group)) +
	geom_abline(intercept=0,slope=1,size=2,alpha=1) +
	#coord_cartesian(xlim=c(5,25),ylim=c(5,25)) +
	scale_color_manual(values=c("#228B22", "#8B4513", "#e60000", "#b8b8b8")) + 
	theme_classic() + facet_wrap(~models,scale="free")
LONG <- PRISMA %>% pivot_longer(cols = 4:9, names_to="models", values_to="brainage")
LONG <- LONG %>% pivot_longer(cols = grep("_CLUSTER",names(LONG)), names_to="developmental", values_to="group")  
LONG$developmental<-gsub("_CLUSTER","",LONG$developmental) ; LONG<-LONG[which(LONG$models==LONG$developmental),] 
LONG$models<-factor(LONG$models, levels = c("DBN","rDBN","tDBN","GTB","rGTB","tGTB")) 
ggplot() +
	geom_smooth(data=LONG,aes(x=age,y=brainage,group=sub,color=group),method = "lm",se=F,size=1.5,alpha=1) +
	geom_point(data=LONG,aes(x=age,y=brainage,color=group),size=1.5,alpha=0.75) +
	geom_abline(intercept=0,slope=1,size=2,alpha=1) +
	#coord_cartesian(xlim=c(5,25),ylim=c(5,25)) +
	scale_color_manual(values=c("#228B22", "#e60000", "#b8b8b8", "#8B4513")) + 
	theme_classic() + facet_wrap(~models,scale="free")

##### Visualize Developmental Trajectories that Parcelated Groups 
PRISMA<-PRISMA[!duplicated(PRISMA$sub),]
OUT<-as.data.frame(names(PRISMA)[4:9]) ; names(OUT)<-"models"
LONG <-pivot_longer(PRISMA,cols=grep("_GROUP",names(PRISMA))[1:6],names_to="models",values_to="braingrps") 
LONG$braingrps<-factor(LONG$braingrps,levels=c("Delayed","Normal","Accelerated"))
LONG$braingrp<-factor(LONG$braingrps,levels=c("Delayed","Accelerated"))
LONG_NEST <- LONG %>% group_by(models) %>% nest # 48:52 / 53:57 
for (COL in c(3,grep("_TIV",names(PRISMA)))){
	MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lm(PRISMA[,COL] ~ braingrp + age + sex + EulerNumber, data = .)))
	MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))  
	for (MOD in 1:dim(MODELS)[1]){
		print(paste0(MODELS$models[MOD],"x",names(PRISMA)[COL])) ; print(model_parameters(MODELS$fit1[[MOD]]))
		OUT[which(OUT$models==gsub("_GROUP","",MODELS$models[MOD])),names(PRISMA)[COL]]<-round(model_parameters(MODELS$fit1[[MOD]])[4,6],digits=3) 
	}
}

### Prepare the Output DataFrame for Visuals
for (PRED in c(53:57)){
	for (MODL in c(4:11)){
		GROUPVAR<-which(names(PDS)==paste0(names(PDS)[MODL],"_GROUP"))
		CONTENT<-cbind(names(PDS)[MODL],names(PDS)[PRED],summarySE(PDS,measurevar=names(PDS)[PRED],groupvars=names(PDS)[GROUPVAR],na.rm=TRUE))
		names(CONTENT)[1:5]<-c("model","predictor","group","N","score") ; CONTENT$score<-scale(CONTENT$score,scale=FALSE)
		if (PRED==53 && MODL==4){
			OUTPUT<-CONTENT
		}
		OUTPUT<-rbind(OUTPUT,CONTENT)
	}
}

### Visualize Differences Between Developmental Groups by Each Predictor
OUTPUT$group<-factor(OUTPUT$group,levels=c("Delayed","Accelerated"))
OUTPUT$model<-factor(OUTPUT$model,levels=c("DBN","rDBN","tDBN","GTB","rGTB","tGTB","hGTBgam","hGTBlmer"))
OUTPUT$predictor<-factor(gsub("_TIV","",OUTPUT$predictor),levels=c("ListSor","PatComp","PictSeq","PictVoc","Flanker"))
pd <- position_dodge(0.50) ; OUTPUT<-OUTPUT[!is.na(OUTPUT$group),]
ggplot(OUTPUT, aes(x=predictor,y=score,colour=group)) + 
	geom_errorbar(aes(ymin=score-se, ymax=score+se), position=pd) +
	geom_point(position=pd,aes(size=0.5)) + geom_line(position=pd) +
	theme_classic() + facet_wrap(~model, scales="free") + 
	#scale_color_manual(values=c("#228B22", "#8B4513", "#e60000"))  
	scale_color_manual(values=c("#228B22", "#e60000"))

######
### Analyze Within-Subject Changes via Person-Specific Centering
######

#Center Each Person By their Specific Mean
WITHIN<-PDS
for (SUB in unique(WITHIN$sub)){
	if (length(which(WITHIN$sub == SUB)) == 1){
		WITHIN<-WITHIN[-c(which(WITHIN$sub == SUB)),]
	} else {
		for (COL in c(3:11,30)){
			WITHIN[which(WITHIN$sub == SUB),COL]<-(WITHIN[which(WITHIN$sub == SUB),COL]-mean(WITHIN[which(WITHIN$sub == SUB),COL]))
		}
	}
}

#Run the Models for Each Brain Age Prediction
LONG <- WITHIN %>% pivot_longer(cols = 4:11, names_to="models", values_to="brainages") ; LONG_NEST <- LONG %>% group_by(models) %>% nest()
MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lmer(brainages ~ age + sex + scanner + EulerNumber + INR + (age | sub), data = ., control=lmerControl(optimizer="optimx",optCtrl=list(method='nlminb')))))
MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))
for (MOD in 1:dim(MODELS)[1]){
	print(MODELS$models[MOD])
	print(model_parameters(MODELS$fit1[[MOD]]))
}

#Visuals of Within-Subject Relationships
FIXED <- ggpredict(MODELS$fit1[[1]], terms = "INR")
RANDOM <- ggpredict(MODELS$fit1[[1]], terms = c("INR", "sub"), type = "random")
ggplot() +
	geom_line(mapping=aes(x=x, y=predicted, group=group),data=RANDOM, alpha=.25, color="#000000") +
	geom_line(mapping=aes(x=x, y=predicted, group=group),data=FIXED, size=3.5, color="#0022ff") +
	theme_classic()

FIXED <- ggpredict(MODELS$fit1[[2]], terms = "INR")
RANDOM <- ggpredict(MODELS$fit1[[2]], terms = c("INR", "sub"), type = "random")
ggplot() +
	geom_line(mapping=aes(x=x, y=predicted, group=group),data=RANDOM, alpha=.25, color="#000000") +
	geom_line(mapping=aes(x=x, y=predicted, group=group),data=FIXED, size=3.5, color="#00d9ff") +
	theme_classic()

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######