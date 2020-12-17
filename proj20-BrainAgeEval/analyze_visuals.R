#!/usr/bin/env Rscript
######################

library(mgcv)
library(lme4)
library(ggplot2)
library(cowplot)
library(corrplot)
library(RColorBrewer)

#########
### Prepare Data By Relabeling Identifiers
#########

data<-read.csv("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/n834_DataFreeze_20201213.csv")
GetTimepoint <- function(TPnum){
	SUBJECTS<-unique(data$sub)
	newdata=data.frame()
	for (x in SUBJECTS){	
		row<-which(data$sub==x)[TPnum]
		addrow<-data[row,]
		newdata<-rbind(newdata,addrow)
		newdata<-newdata[complete.cases(newdata$sub),]
	}
	return(newdata)
}
TP1<-GetTimepoint(1) ; TP1$Session<-1 #213
TP2<-GetTimepoint(2) ; TP2$Session<-2 #204
TP3<-GetTimepoint(3) ; TP3$Session<-3 #189
TP4<-GetTimepoint(4) ; TP4$Session<-4 #143
TP5<-GetTimepoint(5) ; TP5$Session<-5 #84
data<-rbind(TP1,TP2,TP3,TP4,TP5) ; data<-data[order(data$age),]
data$Subject<- 0 ; data$sex<-as.factor(data$sex)
for (x in 1:length(unique(data$sub))){
	subid<-unique(data$sub)[x]
	data[which(data$sub==subid),dim(data)[2]]<-x
}
data$Subject<-as.factor(data$Subject) ; data$Session<-as.factor(data$Session)
data$DBN_MAE<-(data$brainage_DBN - data$age)
data$GTB_MAE<-(data$brainage_GTB - data$age)

#########
### F1: Sample Distribution * Age-related Changes of Covariates
#########

#Age Distribution Ordered by Scan Date
ggplot(data,aes(x=age,y=Subject,group=Subject,color=sex)) + geom_line(size=.8) + geom_point() + scale_color_manual(values=c("#2d81f7","#e62929")) + theme_classic() 
ggplot(data,aes(x=age,y=Subject,group=Subject,color=Session)) + geom_line(size=.8) + geom_point() + scale_color_manual(values=c("#581845","#900C3F","#C70039","#FF5733","#FFC300")) + theme_classic()

#Correlation Matrix of All Numeric Predictors
mat<-data[,c("age","brainage_DBN","brainage_GTB","EulerNumber","T1Income_to_Need","MDDCore","INTLdim","EXTLdim","TtlPrb_Z_sum")]
#complete<-mat[which(complete.cases(mat)==TRUE),] 
Matrix<-cor(mat, use="pairwise.complete.obs")
corrplot.mixed(Matrix, lower.col = "black", number.cex = 1.75)

#Spaghetti Plots of Each Numeric Covariate
ggplot(data, aes(x=age,y=T1Income_to_Need,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic()  
ggsave("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/Figures/Spaghetti_T1Income_to_Need.pdf",plot = last_plot(),device = NULL)

ggplot(data, aes(x=age,y=MDDCore,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic()  
ggsave("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/Figures/Spaghetti_MDD.pdf",plot = last_plot(),device = NULL)

ggplot(data, aes(x=age,y=INTLdim,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic()  
ggsave("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/Figures/Spaghetti_INTLdim.pdf",plot = last_plot(),device = NULL)

ggplot(data, aes(x=age,y=EXTLdim,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic()  
ggsave("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/Figures/Spaghetti_EXTLdim.pdf",plot = last_plot(),device = NULL)

#########
### Subject-Level Trajectories of Brain Age
#########

DBN<-data[,c("sub","ses","age","sex","MDDCore","EulerNumber","brainage_DBN","DBN_MAE")] 
DBN$Model<-"DBN Model" ; names(DBN)[7]<-"brainage" ; names(DBN)[8]<-"MAE"
GTB<-data[,c("sub","ses","age","sex","MDDCore","EulerNumber","brainage_GTB","GTB_MAE")] 
GTB$Model<-"GTB Model" ; names(GTB)[7]<-"brainage" ; names(GTB)[8]<-"MAE"
stack<-rbind(DBN,GTB) ; stack$ses<-as.factor(stack$ses)

### Complete Sample

ggplot(stack, aes(x=age,y=brainage,group=sub,color=Model)) + stat_smooth(method = "lm",se=F,size=0.6,alpha=0.5) + xlim(0,55) + ylim(0,55) + geom_abline(intercept=0,slope=1) + scale_color_manual(values=c("#0022ff","#7c3f00")) + theme_classic() + facet_wrap(~Model)

### Split Scatterplots By Scanner
stack$ses<-as.numeric(stack$ses) ; TRIO<-stack[ which(stack$ses < 4),] ; PRISMA<-stack[ which(stack$ses > 3),] ; stack$ses<-as.numeric(stack$ses)
ggplot(TRIO, aes(x=age,y=brainage,group=sub,color=Model)) + stat_smooth(method = "lm",se=F,size=0.6,alpha=0.5) + xlim(0,55) + ylim(0,55) + geom_abline(intercept=0,slope=1) + scale_color_manual(values=c("#0022ff","#7c3f00")) + theme_classic() + facet_wrap(~Model)
ggplot(PRISMA, aes(x=age,y=brainage,group=sub,color=Model)) + stat_smooth(method = "lm",se=F,size=0.6,alpha=0.5) + xlim(0,55) + ylim(0,55) + geom_abline(intercept=0,slope=1) + scale_color_manual(values=c("#0022ff","#7c3f00")) + theme_classic() + facet_wrap(~Model)

### Split Boxplots for Each Scanner
stack$ses<-as.factor(data$ses) ; stack$sub<-as.factor(data$sub)
ggplot(stack, aes(x=ses, y=MAE,color=Model)) + geom_boxplot() + theme_classic() + scale_color_manual(values=c("#0022ff","#7c3f00")) + theme_classic() + facet_wrap(~Model)

#########
### Run Mixed Effect Models
#########

REFINED<-data[complete.cases(data[,c("sub","ses","MDDCore","EulerNumber","age","sex","brainage_GTB","brainage_GTB")]),]
model_DBN<-lme(brainage_DBN ~ age + sex + MDDCore + EulerNumber, random = ~ 1 | sub, data=REFINED) ; summary(model_DBN)
model_GTB<-lme(brainage_GTB ~ age + sex + MDDCore + EulerNumber, random = ~ 1 | sub, data=REFINED) ; summary(model_GTB)

#########
### Examine Confounding Effects of Image Quality using Euler Number
#########

sensitivity_qa<-data[which(data$EulerNumber > -169),]
DBN<-sensitivity_qa[,c("sub","ses","age","brainage_DBN","DBN_MAE")] 
DBN$Model<-"DBN Model" ; names(DBN)[4]<-"brainage" ; names(DBN)[5]<-"MAE"
GTB<-sensitivity_qa[,c("sub","ses","age","brainage_GTB","GTB_MAE")] 
GTB$Model<-"GTB Model" ; names(GTB)[4]<-"brainage" ; names(GTB)[5]<-"MAE" 
refined_stack<-rbind(DBN,GTB) ; refined_stack$ses<-as.factor(refined_stack$ses)

### Visualize of Image Quality
ggplot(data, aes(x=EulerNumber))+geom_histogram(aes(y=..density..,bins=30), colour="black", fill="white") + geom_density(alpha=.6, fill="#999493")+theme_classic()

### Exclude Bad Quality Scans and Create Boxplots for Each Wave
ggplot(refined_stack, aes(x=ses, y=MAE,color=Model)) + geom_boxplot() + theme_classic() + scale_color_manual(values=c("#0022ff","#7c3f00")) + theme_classic() + facet_wrap(~Model)

### Create Overall Figure with Different Color Scheme
data$AgeGap_DBN<-(data$brainage_DBN - data$age) ; data$AgeGap_GTB<-(data$brainage_GTB - data$age)
data$Euler_INCLUSION<-0 ; data[which(data$EulerNumber > -169),"Euler_INCLUSION"]<-1 
data$Euler_INCLUSION<-as.factor(data$Euler_INCLUSION)
ggplot(data,aes(x=ses,y=AgeGap_GTB,group=sub,color=Euler_INCLUSION)) + geom_point() + scale_color_manual(values=c("#f20000","#000000")) + stat_smooth(method = "lm",se=F,size=0.2,alpha=0.2) + geom_abline(intercept=0,slope=0) + theme_classic()
ggplot(data,aes(x=ses,y=AgeGap_DBN,group=sub,color=Euler_INCLUSION)) + geom_point() + scale_color_manual(values=c("#f20000","#000000")) + stat_smooth(method = "lm",se=F,size=0.2,alpha=0.2) + geom_abline(intercept=0,slope=0) + theme_classic()

### Create Overall Figure with Different Color Scheme

ggplot() + 
	geom_point(data,mapping=aes(x=ses,y=AgeGap_GTB,colour="#7c3f00"),size=.9)+
	geom_point(data,mapping=aes(x=ses,y=AgeGap_DBN,colour="#0022ff"),size=.9)+
	stat_smooth(data,mapping=aes(x=ses,y=AgeGap_GTB,colour="#7c3f00"),method = "lm",se=F,size=.7) +
	stat_smooth(data,mapping=aes(x=ses,y=AgeGap_DBN,colour="#0022ff"),method = "lm",se=F,size=.7) +
	scale_color_manual(values=c("#0022ff","#7c3f00"))+ 
	geom_abline(intercept=0,slope=0,,size=.3,alpha=.3) +
	theme_classic() + facet_wrap(~sub) 

### Create Overall Figure with Different Color Scheme

ggplot() +
	geom_density(TP4,mapping=aes(x=MDDCore,y=..density..), fill="#0eb011",alpha=0.5) +
	geom_density(TP1,mapping=aes(x=MDDCore,y=..density..), fill="#ff0000",alpha=0.5) +
	geom_density(TP2,mapping=aes(x=MDDCore,y=..density..), fill="#0004ff",alpha=0.5) +
	geom_density(TP3,mapping=aes(x=MDDCore,y=..density..), fill="#00d5e0",alpha=0.5) +
	geom_density(TP5,mapping=aes(x=MDDCore,y=..density..), fill="#e6de00",alpha=0.5) +
	theme_classic()

### Create Overall Figure with Different Color Scheme

ggplot(data, aes(x=age,y=MDDCore,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic() 

data$MDD_slope<-NA ; data$MDD_intercept<-NA
for (SUB in  unique(data$sub)){
	temp<-data[which(data$sub == SUB),]
	m1<-lm(temp$age~temp$MDDCore)
	data[which(data$sub==SUB),"MDD_intercept"]<-m1$coefficients[1]
	data[which(data$sub==SUB),"MDD_slope"]<-m1$coefficients[2]
}

ggplot(data, aes(x=MDD_slope))+geom_histogram(aes(y=..density..,bins=30), colour="black", fill="white") + geom_density(alpha=.6, fill="#999493")+theme_classic()


DBN<-data[,c("sub","ses","age","sex","MDDCore","EulerNumber","brainage_DBN","DBN_MAE","MDD_slope")] 
DBN$Model<-"DBN Model" ; names(DBN)[7]<-"brainage" ; names(DBN)[8]<-"MAE"
GTB<-data[,c("sub","ses","age","sex","MDDCore","EulerNumber","brainage_GTB","GTB_MAE","MDD_slope")] 
GTB$Model<-"GTB Model" ; names(GTB)[7]<-"brainage" ; names(GTB)[8]<-"MAE"
stack<-rbind(DBN,GTB) ; stack$ses<-as.factor(stack$ses)

ggplot(data,aes(x=ses,y=AgeGap_DBN,group=sub,color=MDD_slope)) + geom_point(size=0.4) + stat_smooth(method = "lm",se=F,size=0.4,alpha=0.8) + geom_abline(intercept=0,slope=0) + scale_colour_gradient2(low = "#940000",mid = "#787878",high = "#060b9e",midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour") + theme_classic()

ggplot(data,aes(x=ses,y=AgeGap_GTB,group=sub,color=MDD_slope)) + geom_point(size=0.4) + stat_smooth(method = "lm",se=F,size=0.4,alpha=0.8) + geom_abline(intercept=0,slope=0) + scale_colour_gradient2(low = "#940000",mid = "#787878",high = "#060b9e",midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour") + theme_classic()

#########
### Examine Confounding Effects of Image Quality using Euler Number
#########

ggplot(TRIO, aes(x=EulerNumber))+geom_histogram(aes(y=..density..,bins=30), colour="black", fill="white") + geom_density(alpha=.6, fill="#999493")+theme_classic()
ggplot(PRISMA, aes(x=EulerNumber))+geom_histogram(aes(y=..density..,bins=30), colour="black", fill="white") + geom_density(alpha=.6, fill="#999493")+theme_classic()

TRIO$Euler_INCLUSION<-0 ; PRISMA$Euler_INCLUSION<-0
TRIO[which(TRIO$EulerNumber > -206),"Euler_INCLUSION"]<-1 
PRISMA[which(PRISMA$EulerNumber > -48.5),"Euler_INCLUSION"]<-1
EXCLUDE<-rbind(TRIO,PRISMA) ; EXCLUDE$Euler_INCLUSION<-as.factor(EXCLUDE$Euler_INCLUSION)

ggplot(EXCLUDE,aes(x=ses,y=AgeGap_GTB,group=sub,color=Euler_INCLUSION)) + geom_point() + scale_color_manual(values=c("#f20000","#000000")) + stat_smooth(method = "lm",se=F,size=0.2,alpha=0.2) + geom_abline(intercept=0,slope=0) + theme_classic()
ggplot(EXCLUDE,aes(x=ses,y=AgeGap_DBN,group=sub,color=Euler_INCLUSION)) + geom_point() + scale_color_manual(values=c("#f20000","#000000")) + stat_smooth(method = "lm",se=F,size=0.2,alpha=0.2) + geom_abline(intercept=0,slope=0) + theme_classic()

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
