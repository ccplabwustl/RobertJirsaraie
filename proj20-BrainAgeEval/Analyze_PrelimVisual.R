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
### Analyze ANOVAs Using Groups of Brain Age from the PDS Dataset
######

###Histogram of Age Distributions
DATAxAGE<-ddply(COMBINE,"DATA", summarize,age=mean(age))
ggplot(COMBINE, aes(x=age, color=DATA, fill=DATA)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=DATAxAGE, aes(xintercept=age, color=DATA),linetype="longdash",size=3.5) + 
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic()

###Longitudinal Distribution Line Plot
PDS<-PDS[order(PDS$age),] ; PDS$Subject<- 0 
for (x in 1:length(unique(PDS$sub))){
	subid<-unique(PDS$sub)[x]
	PDS[which(PDS$sub==subid),dim(PDS)[2]]<-x
} ; PDS<-PDS[order(PDS$sub),]
PDS$DATA<-factor(PDS$DATA,levels=c("PDS_TRIO","PDS_PRISMA"))
PDS[,c("Subject","ses")]<-lapply(PDS[,c("Subject","ses")],factor)
ggplot() + 
	geom_line(PDS,mapping=aes(x=age, y=Subject, group=Subject, color=DATA),size=1.325,alpha=0.75) +
	geom_point(PDS,mapping=aes(shape=DATA,x=age, y=Subject,color=DATA,fill=DATA),size=3.5) + 
	scale_shape_manual(values=c(16,17))+
	scale_color_manual(values=c("#007A39","#0080FE"))+
	scale_fill_manual(values=c("#007A39","#0080FE"))+
	theme_classic() 

###Histogram of Cogntiive Scores
TEMP <-pivot_longer(COMBINE,cols=16:25,names_to="DOMAINS",values_to="SCORES") 
ggplot(TEMP, aes(x=SCORES, color=DATA, fill=DATA)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.50) +
	scale_color_manual(values=c("#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~DOMAINS, ncol = 5)

###Scatterplot of Age-Related Changes in Cognition
ggplot() +
	geom_point(data=TEMP,aes(x=age,y=SCORES,color=DATA),size=1.5,alpha=0.25) +
	geom_smooth(data=TEMP,aes(x=age,y=SCORES,group=DATA,color=DATA),method = "lm",size=2.4,fullrange=T,geom='line',alpha=0.5,se=FALSE) +
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~DOMAINS, ncol = 5)

###Look for Significant Group Differences by Dataset
LONG <- pivot_longer(COMBINE, cols = 16:25, names_to="models", values_to="brainages") ; LONG_NEST <- LONG %>% group_by(models) %>% nest()
MODELS <- LONG_NEST %>% mutate(fit1 = map(data, ~ lm(brainages ~ age, data = .))) 
MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy)) #MODELS$tidy
for (MOD in 1:dim(MODELS)[1]){
	print(MODELS$models[[MOD]])
	print(paste0(model_parameters(MODELS$fit1[[MOD]])[2,6],"x",model_parameters(MODELS$fit1[[MOD]])[2,8]))
}

###Visualize the Age-Dependency Artifact for all Models
model_parameters(gamm4(BRAINAGE_SD_NORM ~ s(age) + sex + EulerNumber, random=~(1 | sub), data=PRISMA))
model_parameters(gamm4(BRAINAGE_SD_NORM ~ s(age) + sex + EulerNumber, random=~(1 | sub), data=TRIO))
model_parameters(gam(BRAINAGE_SD_NORM ~ s(age) + sex + EulerNumber + numNavs_T1w, data=HCP))  
LONG <-pivot_longer(COMBINE,cols=grep("_aERROR",names(COMBINE))[1:6],names_to="MODELS",values_to="BAGS")
ggplot(LONG, aes(x=age,y=BAGS,group=MODELS,color=MODELS)) + 
	geom_point(size=1.5,alpha=0.5) +
	theme_classic() + facet_wrap(~DATA, scale="free") + 
	stat_smooth(method = "lm",se=F,size=2.5,alpha=1,fullrange=T) +
	scale_color_manual(values=c("#0022ff","#00d9ff","#06029e","#7c3f00","#db9e65","#472309"))



