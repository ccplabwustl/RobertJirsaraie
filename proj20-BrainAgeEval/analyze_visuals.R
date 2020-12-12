#!/usr/bin/env Rscript
######################

library(mgcv)
library(ggplot2)
library(cowplot)
library(corrplot)
library(RColorBrewer)

#########
### Prepare Data By Relabeling Identifiers
#########

data<-read.csv("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/n834_DataFreeze_20201122.csv")
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

#########
### Distribution of Scan Age
#########

ggplot(data,aes(x=age,y=Subject,group=Subject,color=sex)) + geom_line(size=.8) + geom_point() + scale_color_manual(values=c("#2d81f7","#e62929")) + theme_classic() 
ggplot(data,aes(x=age,y=Subject,group=Subject,color=Session)) + geom_line(size=.8) + geom_point() + scale_color_manual(values=c("#581845","#900C3F","#C70039","#FF5733","#FFC300")) + theme_classic()
ggplot(data,aes(x=age,y=Subject,group=Subject,color=Session)) + geom_line(size=.8) + geom_point() + scale_color_manual(values=c("#173F5F","#20639B","#3CAEA3","#F6D55C","#ED553B")) + theme_classic() 

#########
### Correlation Matrix of All Numeric Predictors
#########

mat<-data[,c("age","brainage_DBN","T1Income_to_Need","MDDCore","INTLdim","EXTLdim","TtlPrb_Z_sum")]
complete<-mat[which(complete.cases(mat)==TRUE),] 
Matrix<-cor(complete, use="pairwise.complete.obs")
corrplot.mixed(Matrix, lower.col = "black", number.cex = 1.75)

#########
### Spaghetti Plots of Each Numeric Covariate
#########

ggplot(data, aes(x=age,y=T1Income_to_Need,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic()  
ggsave("/scratch/rjirsara/proj20-BrainAgeEval/Analysis/Figures/Spaghetti_T1Income_to_Need.pdf",plot = last_plot(),device = NULL)

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

ggplot(data, aes(x=age,y=brainage_DBN,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + xlim(0, 22) + ylim(0,50) + geom_abline(intercept=0,slope=1) + theme_classic()
ggplot(data, aes(x=age,y=brainage_DBN,group=sub)) + geom_line(col="blue") + geom_point(size=1,alpha=0.8) + xlim(0, 60) + ylim(0,60) + geom_abline(intercept=0,slope=1) + theme_classic()

data$COHORT<-0 ; data[which(data$ses<=3),"COHORT"]<-1 ; data[which(data$ses>=3),"COHORT"]<-2
ggplot(data, aes(x=age,y=brainage_DBN,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + xlim(0, 22) + ylim(0,50) + geom_abline(intercept=0,slope=1) + theme_classic()+facet_wrap(~COHORT)
ggplot(data, aes(x=age,y=brainage_DBN,group=sub)) + geom_line(col="blue") + geom_point(size=1,alpha=0.8) + xlim(0, 60) + ylim(0,60) + geom_abline(intercept=0,slope=1) + theme_classic()+facet_wrap(~COHORT)

#########
### Within-Subject Effects on Mean Absolute Error and Psychopathology
#########

data$MDDCorez<-scale(data$MDDCore)
data$brainage_DBNz<-scale(data$brainage_DBN)
data$brainage_MAE_DBN<-abs(data$age-data$brainage_DBN)
data$ses<-as.factor(data$ses) ; data$sub<-as.factor(data$sub)
ggplot(data, aes(x=ses, y=brainage_MAE_DBN)) + geom_boxplot() + theme_classic()
TP1<-data[which(data$ses==1),] ; summary(TP1$brainage_MAE_DBN) ; sd(TP1$brainage_MAE_DBN) # 3.52-2.8
TP2<-data[which(data$ses==2),] ; summary(TP2$brainage_MAE_DBN) ; sd(TP2$brainage_MAE_DBN) # 4.11-3.3
TP3<-data[which(data$ses==3),] ; summary(TP3$brainage_MAE_DBN) ; sd(TP3$brainage_MAE_DBN) # 4.82-4.1
TP4<-data[which(data$ses==4),] ; summary(TP4$brainage_MAE_DBN) ; sd(TP4$brainage_MAE_DBN) # 22.52-5.33
TP5<-data[which(data$ses==5),] ; summary(TP5$brainage_MAE_DBN) ; sd(TP5$brainage_MAE_DBN) # 23.24-5.05

####

ggplot(data, aes(x=age,y=brainage_MAE_DBN,group=sub)) + geom_point() + geom_point(group="sub") + stat_smooth(method = "lm",se = FALSE) + facet_wrap(~sub)
gamm4(formula = brainage ~ age + sex + incometoneeds + MDD, random=as.formula(~(age|sub)), data=dataSubj, REML=T)$gam
model_CorticalGray <-lmer(ZCortGrayVol~ 1+Zage_cntr13+Zage_cntr13_sqr+Zsex+ZIntraCranialVol_cntr+Zincome_cntr+Zage_cntr13*Zincome_cntr + (Zage_cntr13 | Subid), data=pdslong_stan) summary(model_CorticalGray)
r2beta(model_CorticalGray, partial = TRUE, method = "nsj", data = NULL)
~(age|sub) #differenece as a function of age
~(sub|age)
ggplot(data, aes(x=age,y=TtlPrb_Z_sum,group=sub)) + stat_smooth(method = "lm",se=F,size=0.5,alpha=0.5) + theme_classic()  
ggsave("/scratch/rjirsara/projects/BrainAgeEval/Analysis/Figures/Spaghetti_TtlPrb_Z_sum.pdf",plot = last_plot(),device = NULL)
gamm4(formula = x, random=as.formula(randomFormula), data=dataSubj, REML=T)$gam

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
