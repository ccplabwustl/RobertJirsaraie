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
COMBINE<-read.csv(paste0(NEURO_PATH,"/COMBINE.csv"))
PRISMA<-read.csv(paste0(NEURO_PATH,"/PRISMA.csv"))
TRIO<-read.csv(paste0(NEURO_PATH,"/TRIO.csv"))
HCP<-read.csv(paste0(NEURO_PATH,"/HCP.csv"))

######
### Model Accuracy
######

###Scatterplots of Brain Age Models Fit with Age By Dataset
TEMP <-pivot_longer(COMBINE,cols=4:9,names_to="MODELS",values_to="BRAINAGE") 
TEMP$MODELS<-factor(TEMP$MODELS,levels=c("DBN","rDBN","tDBN","GTB","rGTB","tGTB"))
ggplot(TEMP, aes(x=age,y=BRAINAGE,group=DATA,color=DATA,shape=DATA)) + 
	geom_abline(intercept=0,slope=1,size=2,alpha=1,linetype="dashed") +
	geom_point(size=1.45,alpha=0.575) +
	stat_smooth(method = "lm",se=F,size=2.5,alpha=0.5,fullrange=T) +
	scale_shape_manual(values=c(16,17,15))+
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS,scales = "free")

###Benchmark Developmental Brain Age Metrics
for (COL in 4:11){
	DATA<-PRISMA ; TYPE<-"LONG" ; print(paste0(names(DATA)[COL]))
	if (TYPE == "CROSS"){
		COEF<-round(lm(DATA[,COL]~DATA[,"age"])$coefficients,digits=2)
	} else if (TYPE == "LONG"){
		COEF<-round(fixef(lmer(DATA[,COL] ~ DATA[,"age"] + (1 | DATA[,"sub"]))),digits=2)
	}
	MAE<-round(mean(abs((DATA[,COL]-DATA[,"age"])),na.rm=T),digits=2)
	CORR<-round(cor(DATA[,COL],DATA[,"age"],use = "complete.obs"),digits=2)
	print(paste0("MAEx",MAE,"_SLOPEx",COEF[2],"_Interceptx",COEF[1],"_CORRx",CORR))
}

######
### Model Sensitivity
######

###Individual Differences in Cognition #3:8 / 23:28 / 29:34 
OUT<-as.data.frame(names(HCP)[3:8]) ; names(OUT)[1]<-"models" 
TEMP <-pivot_longer(HCP,cols=3:8,names_to="models",values_to="brainages") 
TEMP_NEST <- TEMP %>% group_by(models) %>% nest
for (COL in grep("_RAW",names(HCP))[1:5]){
	MODELS <- TEMP_NEST %>% mutate(fit1 = map(data, ~ lm(HCP[,COL] ~ brainages + age + sex + EulerNumber, data = .)))
	MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))  
	for (MOD in 1:dim(MODELS)[1]){
		print(paste0(MODELS$models[MOD],"x",names(HCP)[COL])) 
		print(model_parameters(MODELS$fit1[[MOD]])) #print(lm.beta(MODELS$fit1[[MOD]]))
		OUT[which(OUT$models==MODELS$models[MOD]),names(HCP)[COL]]<-round(model_parameters(MODELS$fit1[[MOD]])[2,6],digits=3) 
	}
}

OUT<-as.data.frame(names(PRISMA)[4:9]) ; names(OUT)[1]<-"models" #4:9 / 32:37 / 38:43 
TEMP <-pivot_longer(PRISMA,cols=4:9,names_to="models",values_to="brainages") 
TEMP_NEST <- TEMP %>% group_by(models) %>% nest
for (COL in grep("_RAW",names(PRISMA))[1:5]){
	MODELS <- TEMP_NEST %>% mutate(fit1 = map(data, ~ lmer(PRISMA[,COL] ~ brainages + age + sex + EulerNumber + (1 | sub), data = .)))
	MODELS <- MODELS %>% mutate(tidy = map(fit1, broom::tidy))  
	for (MOD in 1:dim(MODELS)[1]){
		print(paste0(MODELS$models[MOD],"x",names(PRISMA)[COL])) 
		print(model_parameters(MODELS$fit1[[MOD]])) ; #print(lm.beta.lmer(MODELS$fit1[[MOD]]))
		OUT[which(OUT$models==MODELS$models[MOD]),names(PRISMA)[COL]]<-round(model_parameters(MODELS$fit1[[MOD]])[2,6],digits=3) 
	}
}

######
### Create Scatterplot to Show the Relability of Brain Age Perdictions at an Subject Level
######

ROWS1 <- as.numeric(row.names(COMBINE[order(COMBINE[which(COMBINE$DATA=="PDS_TRIO"),"age"]),]))
ROWS2 <- as.numeric(row.names(COMBINE[order(COMBINE[which(COMBINE$DATA=="PDS_PRISMA"),"age"]),]))+dim(TRIO)[1]
ROWS3 <- as.numeric(row.names(COMBINE[order(COMBINE[which(COMBINE$DATA=="HCP_TEST"),"age"]),]))+(dim(TRIO)[1]+dim(PRISMA)[1])
COMBINE<-COMBINE[c(ROWS1,ROWS2,ROWS3),] ; COMBINE[,"ORDER"]<-1:nrow(COMBINE)
#Scatterplots of MAEs
TEMP <-pivot_longer(COMBINE,cols=38:43,names_to="MODELS",values_to="BRAINAGE") #4:9 - 32:37 - 38:43
TEMP$MODELS<-factor(gsub("_NORM","",gsub("_ERROR","", TEMP$MODELS)),levels=c("DBN","rDBN","tDBN","GTB","rGTB","tGTB"))
ggplot(TEMP, aes(x=ORDER, y=BRAINAGE, color=MODELS)) + 
	geom_rect(aes(xmin=0,xmax=432,ymin=-Inf,ymax=Inf),alpha=1,fill="#70a078") +
	geom_rect(aes(xmin=433,xmax=713,ymin=-Inf,ymax=Inf),alpha=1,fill="#70a6f8") +
	geom_rect(aes(xmin=713,xmax=1106,ymin=-Inf,ymax=Inf),alpha=1,fill="#9d6953") + 
	geom_hline(yintercept=0, color="#000000", size=1.5) + geom_point(size=1.75) + theme_classic() +
	scale_color_manual(values=c("#0022ff","#00d9ff","#040080","#7c3f00","#db9e65","#472309"))

#Scatterplots of MAEs REDUCED
INDEX<-c((28*1:10),(28*1:10)+432,(35*1:10)+713)
TEMP <-pivot_longer(COMBINE,cols=38:43,names_to="MODELS",values_to="BRAINAGE") #4:9 - 32:37 - 38:43
TEMP$MODELS<-factor(gsub("_NORM","",gsub("_ERROR","",TEMP$MODELS)),levels=c("DBN","rDBN","tDBN","GTB","rGTB","tGTB"))
ggplot(TEMP, aes(x=ORDER, y=BRAINAGE, color=MODELS)) + 
	geom_rect(aes(xmin=0,xmax=432,ymin=-Inf,ymax=Inf),alpha=1,fill="#70a078") +
	geom_rect(aes(xmin=433,xmax=713,ymin=-Inf,ymax=Inf),alpha=1,fill="#70a6f8") +
	geom_rect(aes(xmin=713,xmax=1106,ymin=-Inf,ymax=Inf),alpha=1,fill="#9d6953") + 
	geom_hline(yintercept=0, color="#000000", size=2.5) + geom_point(size=1.75) + theme_classic() +
	scale_color_manual(values=c("#0022ff","#00d9ff","#040080","#7c3f00","#db9e65","#472309"))

#Histogram of Reliability by Models
HIST_RAW <-pivot_longer(COMBINE,cols=c(59,61),names_to="MODELS",values_to="DEVS") 
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_SD_RAW), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_SD_RAW), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_RAW, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")
HIST_NORM <-pivot_longer(COMBINE,cols=c(60,62),names_to="MODELS",values_to="DEVS")
MEANS1<-ddply(COMBINE,"DATA", summarize, mean=mean(GTB_SD_NORM), model="GTB") 
MEANS2<-ddply(COMBINE,"DATA", summarize, mean=mean(DBN_SD_NORM), model="DBN") 
MEAN<-rbind(MEANS1,MEANS2)
ggplot(HIST_NORM, aes(x=DEVS, color=MODELS, fill=MODELS)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=MEAN, aes(xintercept=mean, color=model),linetype="longdash",size=1.5) + 
	scale_color_manual(values=c("#0022ff","#0022ff","#7c3f00","#7c3f00")) +
	scale_fill_manual(values=c("#0022ff","#7c3f00")) +
	theme_classic() + facet_wrap(~DATA,scale="free")

#Histogram of Reliability by Test Sample
DATAxERRORxNORM<-ddply(COMBINE,"DATA", summarize, error=mean(GTB_aSD_NORM))
ggplot(COMBINE, aes(x=GTB_aSD_NORM, color=DATA, fill=DATA)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.70) +
	geom_vline(data=DATAxERRORxNORM, aes(xintercept=error, color=DATA),linetype="longdash",size=3.5) + 
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic()

###Get Model Parameters Regarding Individual Differences in Relability
model_parameters(lmer(BRAINAGE_SD_NORM ~ age + sex + EulerNumber + (1 | sub),data=PRISMA))
model_parameters(lmer(BRAINAGE_SD_NORM ~ age + sex + EulerNumber +  (1 | sub),data=TRIO))
model_parameters(lm(BRAINAGE_SD_NORM ~ age + sex + EulerNumber + numNavs_T1w, data=HCP)) 
LONG <-pivot_longer(COMBINE,cols=grep("_SD_NORM",names(COMBINE))[2:3],names_to="MODELS",values_to="DEVS")
ggplot(LONG, aes(x=age,y=DEVS,group=DATA,color=DATA,shape=DATA)) + 
	geom_point(size=1.5,alpha=0.7) + 	scale_shape_manual(values=c(16,17,15)) +
	stat_smooth(method = "loess",se=F,size=3,alpha=1,fullrange=T) +
	scale_color_manual(values=c("#007A39","#0080FE","#A0522D")) +
	scale_fill_manual(values=c("#007A39","#0080FE","#A0522D")) +
	theme_classic() + facet_wrap(~MODELS)

###Create Visuals of Reliability Differences in Age
NEURO_PATH<-"/Users/Jirsaraie/Box/Research/proj20-BrainAgeEval/audits"
LONG<-read.csv(paste0(NEURO_PATH,"/NEURO_Traj.csv"))
FINAL<-read.csv(paste0(NEURO_PATH,"/NEURO_Grp.csv"))
ggplot(LONG, aes(x=age,y=Values, group=ROIs,color=FEATURE)) + 
	geom_smooth(method = "loess",se=F,size=0.2,alpha=1,fullrange=F) + 
	scale_color_manual(values=c("#0062ff","#f70505","#28b03f","#ffa600")) + 
	theme_classic() + facet_wrap(~DATA, scale="free")

#Variability Between Age Groups By Feature
ggplot(FINAL, aes(x=age_bin,y=DEV_MIN_COV,group=feature,color=FEATURE)) +
	geom_line(size=0.6,alpha=0.5) + geom_point(size=0.7,alpha=0.6) +
	geom_point(FINAL,mapping=aes(x=age_bin,y=MEAN_MIN_COV,color="#000000"),size=3,shape=19) +
	scale_color_manual(values=c("#000000","#0062ff","#f70505","#28b03f","#ffa600")) + 
	theme_classic() + facet_wrap(~data,scale="free")

#Variability Between Age Groups By Feature Type
FINAL$age_bin<-factor(FINAL$age_bin,levels=c("YOUNG","MID","OLD"))
ggplot(FINAL, aes(age_bin, DEV_RAW_COV)) + 
    geom_jitter(width=0.4, alpha=1,size=0.5,aes(colour = FEATURE)) + 
    geom_point(FINAL,mapping=aes(x=age_bin, y=MEAN_RAW_COV),size=3,shape=19) + 
    scale_color_manual(values=c("#0062ff","#f70505","#28b03f","#ffa600")) + 
    facet_wrap(~data,scale="free") + theme_classic() + ylim(0.14,0.20) 

#Variability Between Age Groups Stats
anova(lmer(DEV_MIN_COV~age_bin + (1|feature),PDSTrio))
anova(lmer(DEV_MIN_COV~age_bin + (1|feature),PDSPrisma))
anova(lmer(DEV_MIN_COV~age_bin + (1|feature),HCPTest))

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######