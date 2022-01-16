#!/usr/bin/env Rscript
######################

set.seed(777)
library(caret)
library(mgcv)
library(lme4)
library(dplyr)
library(ggplot2)
library(cowplot)
library(cowplot)
library(stringr)
library(corrplot)
library(ggeffects)
library(RColorBrewer)
library(data.table) ; library(tidyverse) ; library(plyr)

PDS<-read.csv("/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/study-PDS/Analysis/n834_DataFreeze_20210121.csv")
HCP<-read.csv("/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/study-HCP/Analysis/n789_DataFreeze_2020211.csv")

######
### Split Data in Half For Retraining Each Model with HCP-D
######

INDEX <- createDataPartition(HCP[,"age"], p=0.5, list=F, times=1) ; INDEX<-INDEX[-c(123),] 
HCP$Split<-NA ; HCP[-INDEX,"Split"]<-0 ; HCP[INDEX,"Split"]<-1 ; HCP$Split<-as.factor(HCP$Split)
HCP<-HCP[sample(nrow(HCP)),] ; row.names(HCP)<-NULL

chisq.test(table(HCP$Split, HCP$sex)) #0.4992
t.test(HCP$EulerNumber~HCP$Split) #0.7857
t.test(HCP$age~HCP$Split) #0.959

HCP$Age<-HCP$age ; HCP$Female<-HCP$sex-1 ; HCP$Euler<-HCP$EulerNumber ; HCP<-HCP[,-c(1:9)] 
write.csv(HCP, "/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/pipes/replicate_devhcp/n789_ReTrainSplit_20210305.csv", row.names=F)

######
### Create Histogram of Euler Numbers and Brain Age Predictions Between Each Dataset
######

tPDS<-PDS[,c("sub","ses","age","sex","EulerNumber","brainage_DBN","brainage_rDBN","brainage_GTB_FSxCross","brainage_rGTB_FSxCross")]
tHCP<-HCP[,1:9] ; tHCP["site"]<-0 ; names(tHCP)[8]<-"ses" ; names(tPDS)[8:9]<-c("brainage_GTB","brainage_rGTB")
tHCP$DATA<-"HCP" ; tPDS$DATA<-"PDS" ; COMBINE<-rbind(tPDS,tHCP) ; COMBINE$DATA<-as.factor(COMBINE$DATA)
COMBINE<-COMBINE[!is.na(COMBINE$EulerNumber),] ; COMBINE$error<-abs(COMBINE$age-COMBINE$brainage_GTB) 
MALES<-COMBINE[COMBINE$sex==1,] ; FEMALES<-COMBINE[COMBINE$sex==2,]
#COMBINE<-COMBINE[!between(COMBINE$ses,1,3),]

mu <- ddply(MALES, "DATA", summarise, grp.mean=mean(error))
ggplot(MALES, aes(x=error, color=DATA, fill=DATA)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.40)+
	geom_vline(data=mu, aes(xintercept=grp.mean, color=DATA),linetype="longdash",size=2)+
	scale_color_manual(values=c("#000000","#b51818", "#000000"))+
	scale_fill_manual(values=c("#000000","#b51818", "#000000"))+
	theme_classic()

######
### Create Scatterplots for Each Facet of Data
######

tCOMBINE<-COMBINE[,c("sub","ses","age","brainage_DBN","brainage_rDBN","brainage_GTB","brainage_rGTB","DATA")]
names(tCOMBINE)<-c("sub","ses","age","DBN","rDBN","GTB","rGTB","DATA") ; tCOMBINE[tCOMBINE$ses==0,"ses"]<-6
ONE<-tPDS[tPDS$ses==1,] ; TWO<-tPDS[tPDS$ses==2,] ; THREE<-tPDS[tPDS$ses==3,] ; TRIO<-tPDS[tPDS$ses<4,]
FOUR<-tPDS[tPDS$ses==4,] ; FIVE<-tPDS[tPDS$ses==5,] ; PRISMA<-tPDS[tPDS$ses>3,] 
tPDS<-tCOMBINE[as.numeric(tCOMBINE$ses)<6,] ; tHCP<-tCOMBINE[as.numeric(tCOMBINE$ses)>5,]

#Cross-sectional Models
LONG <- tCOMBINE %>% pivot_longer(cols = 4:7, names_to = "model", values_to = "brainage")
ggplot(LONG, aes(x=age,y=brainage,group=model,color=model)) + 
	geom_point(size=0.5,alpha=0.25) +
	geom_smooth(method = "lm",se=F,size=1,alpha=0.5,fullrange=TRUE) +
	geom_abline(intercept=0,slope=1,size=1.25,alpha=1) +
	scale_color_manual(values=c("#0022ff","#7c3f00","#00d9ff","#db9e65","#000000")) +
	facet_wrap(~ses) + 
	theme_classic()

#Longitudinal Models
FIXED<-as.data.frame(matrix(0,nrow=5,ncol=1)) ; tPDS<-tCOMBINE[as.numeric(tCOMBINE$ses)<6,]
for (INDEX in 4:7){
	MOD<-lmer(tPDS[,INDEX] ~ age + (1 | sub), data=tPDS) 
	output<-as.data.frame(ggpredict(MOD,"age")[1:2])
	names(output)[2]<-names(tPDS)[INDEX]
	FIXED<-cbind(FIXED,output)
}

FIXED<-FIXED[,-c(1,4,6,8)] ; names(FIXED)[1]<-"age"
LONG <- FIXED %>% pivot_longer(cols = 2:5, names_to = "model", values_to = "brainage")
RANDOM <- tPDS %>% pivot_longer(cols = 4:7, names_to = "model", values_to = "brainage")
ggplot() +  
	geom_point(RANDOM,mapping=aes(x=age,y=brainage,color=model),size=0.85,alpha=0.5) +
	geom_smooth(LONG,mapping=aes(x=age,y=brainage,group=model,color=model),method = "lm",se=F,size=1.75,alpha=0.5,fullrange=TRUE) +
	geom_abline(intercept=0,slope=1,size=2,alpha=1) + xlim(5,23) + ylim(5,60) +
	scale_color_manual(values=c("#0022ff","#7c3f00","#00d9ff","#db9e65","#000000")) +
	theme_classic()

######
### Create Figures to Compare Intercept, Slope, Correlation, and Mean Absolute Errors for each dataset
######

ONE<-tPDS[tPDS$ses==1,] ; TWO<-tPDS[tPDS$ses==2,] ; THREE<-tPDS[tPDS$ses==3,] ; FOUR<-tPDS[tPDS$ses==4,] ; FIVE<-tPDS[tPDS$ses==5,]
TRIO<-tPDS[tPDS$ses<4,] ; PRISMA<-tPDS[tPDS$ses>3,] ; tPDS<-tCOMBINE[as.numeric(tCOMBINE$ses)<6,] ; tHCP<-tCOMBINE[as.numeric(tCOMBINE$ses)>5,]

OUTPUT<-as.data.frame(matrix(0,nrow=1,ncol=4)) ; INDEX<-0
LABELS<-c("PDS-1","PDS-2","PDS-3","PDS-TRIO","PDS-4","PDS-5","PDS-PRISMA","HCP-Total","PDS-Total")
for (DATA in list(ONE, TWO, THREE, TRIO, FOUR, FIVE, PRISMA, tHCP, tPDS)){
	INDEX=((INDEX+1)) ; DATA_LABEL=LABELS[INDEX]
	if (INDEX == 1 || INDEX == 2 || INDEX == 3 || INDEX == 5 || INDEX == 6 || INDEX == 8){
		for (COLNUM in 4:7){
			COEF<-lm(DATA[,COLNUM]~DATA[,"age"])$coefficients
			CORR<-cor(DATA[,COLNUM],DATA[,"age"])
			MAE<-mean(abs((DATA[,COLNUM]-DATA[,"age"])))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"intercept",COEF[1]))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"slope",COEF[2]))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"correlation",CORR))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"error",MAE))
		}
	}
	if (INDEX == 4 || INDEX == 7 || INDEX == 9){
		for (COLNUM in 4:7){
			lmer(DATA[,COLNUM] ~ DATA[,"age"] + (1 | DATA[,"sub"])) 
			COEF<-fixef(lmer(DATA[,COLNUM] ~ DATA[,"age"] + (1 | DATA[,"sub"])))
			CORR<-cor(DATA[,COLNUM],DATA[,"age"])
			MAE<-mean(abs((DATA[,COLNUM]-DATA[,"age"]))) 
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"intercept",COEF[1]))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"slope",COEF[2]))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"correlation",CORR))
			OUTPUT<-rbind(OUTPUT,c(DATA_LABEL,names(DATA)[COLNUM],"error",MAE))
		}
	}
}
names(OUTPUT)<-c("Dataset","Model","Metric","Value") ; OUTPUT<-OUTPUT[-1,] ; OUTPUT$Value<-as.numeric(OUTPUT$Value) ; INDEX<-0
OUTPUT$Dataset<-as.character(OUTPUT$Dataset) ; OUTPUT$Dataset<-factor(OUTPUT$Dataset,levels=unique(OUTPUT$Dataset))
OUTPUT$Value<-round(OUTPUT$Value,digits=3)

SUBSET<-OUTPUT[which(OUTPUT$Metric=="correlation"),] 
ggplot(SUBSET, aes(x=Dataset, y=Value, color=Model, group=Model)) + geom_line(size=1.25) + 
	geom_hline(aes(yintercept = 0)) +
	geom_hline(aes(yintercept = 1)) +
	geom_line() + geom_point(shape=21, fill="white") + 
	geom_text(aes(label=Value),hjust=0, vjust=0) +
	theme_classic() + #ylim(-10,15) + 
	scale_color_manual(values=c("#0022ff","#7c3f00","#00d9ff","#db9e65","#000000"))

######
### Create Overlapping Age Distributions Between HCP-D and PDS And Evaluation Figures
######

COMBINE<-COMBINE[!between(COMBINE$ses,1,3),] 
mu <- ddply(COMBINE, "DATA", summarise, grp.mean=mean(age))
ggplot(COMBINE, aes(x=age, color=DATA, fill=DATA)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.40)+
	geom_vline(data=mu, aes(xintercept=grp.mean, color=DATA),linetype="longdash",size=2)+
	scale_color_manual(values=c("#000000","#b51818", "#000000"))+
	scale_fill_manual(values=c("#000000","#b51818", "#000000"))+
	theme_classic()

COMBINE<-COMBINE[between(COMBINE$age,16,19),]
mu <- ddply(COMBINE, "DATA", summarise, grp.mean=mean(age))
ggplot(COMBINE, aes(x=age, color=DATA, fill=DATA)) +
	geom_histogram(aes(y=..count..), position="identity", alpha=0.40)+
	geom_vline(data=mu, aes(xintercept=grp.mean, color=DATA),linetype="longdash",size=2)+
	scale_color_manual(values=c("#000000","#b51818", "#000000"))+
	scale_fill_manual(values=c("#000000","#b51818", "#000000"))+
	theme_classic()

LONG <- COMBINE %>% pivot_longer(cols = 6:9, names_to = "model", values_to = "brainage")
ggplot(LONG, aes(x=age,y=brainage,group=model,color=model)) + 
	geom_point(size=0.5,alpha=0.25) +
	geom_smooth(method = "lm",se=F,size=1,alpha=0.5,fullrange=TRUE) +
	geom_abline(intercept=0,slope=1,size=1.25,alpha=1) +
	scale_color_manual(values=c("#0022ff","#7c3f00","#00d9ff","#db9e65","#000000")) +
	facet_wrap(~DATA) + 
	theme_classic()

RHCP=COMBINE[COMBINE$DATA=="HCP",] ; RPDS=COMBINE[COMBINE$DATA=="PDS",]
OUTPUT<-as.data.frame(matrix(0,nrow=1,ncol=4))
for (COLNUM in grep("brainage",names(RHCP))){
	FITTED_HCP<-lm(RHCP[,COLNUM]~age,data=RHCP)$coefficients
	FITTED_PDS<-lm(RPDS[,COLNUM]~age,data=RPDS)$coefficients
	OUTPUT<-rbind(OUTPUT,c(unlist(strsplit(names(RHCP)[COLNUM],"_"))[2],"intercept",FITTED_HCP[1],FITTED_PDS[1]))
	OUTPUT<-rbind(OUTPUT,c(unlist(strsplit(names(RHCP)[COLNUM],"_"))[2],"slope",FITTED_HCP[2],FITTED_PDS[2]))
	COR_HCP<-cor(RHCP[,COLNUM],RHCP$age) ; COR_PDS<-cor(RPDS[,COLNUM],RPDS$age)
	OUTPUT<-rbind(OUTPUT,c(unlist(strsplit(names(RHCP)[COLNUM],"_"))[2],"correlation",COR_HCP,COR_PDS))
	MAE_HCP<-mean(abs(RHCP[,COLNUM]-RHCP$age)) ; MAE_PDS<-mean(abs(RPDS[,COLNUM]-RPDS$age))
	OUTPUT<-rbind(OUTPUT,c(unlist(strsplit(names(RHCP)[COLNUM],"_"))[2],"error",MAE_HCP,MAE_PDS))
}
names(OUTPUT)<-c("Model","Metric","HCP","PDS") ; OUTPUT<-OUTPUT[-1,]
OUTPUT[,c(3,4)]<-lapply(OUTPUT[,c(3,4)],as.numeric) 
OUTPUT <- OUTPUT %>% arrange(Metric,Model)

OUTPUT <- OUTPUT %>% pivot_longer(cols = 3:4, names_to = "Dataset", values_to = "Value")
ggplot(OUTPUT, aes(x=Dataset, y=Value, color=Model, group=Model)) + geom_line(size=1.5) + 
    geom_hline(data = subset(OUTPUT, Metric =="intercept"),aes(yintercept = 0)) +
    geom_hline(data = subset(OUTPUT, Metric =="slope"),aes(yintercept = 1)) +
    geom_line() + geom_point(shape=21, fill="white") +  
    facet_wrap(~Metric, scales = "free") +
    theme_classic() +
    scale_color_manual(values=c("#0022ff","#7c3f00","#00d9ff","#db9e65","#000000"))

######
### Tests of Significance Between Model Performance Metrics
######

CHECKMARK

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######