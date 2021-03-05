#!/usr/bin/env Rscript
######################

library(plyr)
library(mgcv)
library(lme4)
library(dplyr)
library(ggplot2)
library(cowplot)
library(corrplot)
library(RColorBrewer)

plot_histogram <- function(df, feature) {
    plt <- ggplot(df, aes(x=eval(parse(text=feature)))) +
    geom_histogram(aes(y = ..density..), alpha=0.7, fill="#33AADE", color="black") +
    geom_density(alpha=0.3, fill="red") +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density") + theme_classic()
    print(plt)
}

plot_multi_histogram <- function(df, feature, label_column){
    plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
    plt + guides(fill=guide_legend(title=label_column)) + theme_classic()
}

data<-read.csv("/scratch/rjirsara/proj20-BrainAgeEval/analysis/n834_DataFreeze_20210101.csv")

#########
### Prepare Slide 1 of Overall Effects
#########

ROIS<-data[,c(4,21:26)] ; colnames(ROIS)<-gsub("brainage_GTB_","",names(ROIS)) ; colnames(ROIS)<-gsub("COMBATx","",names(ROIS))
ROIS<-ROIS[, c("age","FSxCross","GAMxCross","LMERxCross","FSxLong","GAMxLong","LMERxLong")]
corrplot.mixed(cor(ROIS, use="pairwise.complete.obs") , lower.col = "black", number.cex = 1.6)
ROIS<-cbind(data[,1:2],ROIS)

ROIS$MAE_FSxCross<-ROIS$FSxCross-ROIS$age ; ggplot(ROIS, aes(x=MAE_FSxCross))+geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.6, fill="#F05E23")+geom_vline(xintercept = 0, size=2) + xlim(-10, 40) + theme_classic()
ROIS$MAE_GAMxCross<-ROIS$GAMxCross-ROIS$age ; ggplot(ROIS, aes(x=MAE_GAMxCross))+geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.6, fill="#F05E23")+geom_vline(xintercept = 0, size=2) + xlim(-10, 40) + theme_classic()
ROIS$MAE_LMERxCross<-ROIS$LMERxCross-ROIS$age ; ggplot(ROIS, aes(x=MAE_LMERxCross))+geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.6, fill="#F05E23")+geom_vline(xintercept = 0, size=2) + xlim(-10, 40) + theme_classic()

ROIS$MAE_FSxLong<-ROIS$FSxLong-ROIS$age ; ggplot(ROIS, aes(x=MAE_FSxLong))+geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.6, fill="#0b6623")+geom_vline(xintercept = 0, size=2) + xlim(-10, 40) + theme_classic()
ROIS$MAE_GAMxLong<-ROIS$GAMxLong-ROIS$age ; ggplot(ROIS, aes(x=MAE_GAMxLong))+geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.6, fill="#0b6623")+geom_vline(xintercept = 0, size=2) + xlim(-10, 40) + theme_classic()
ROIS$MAE_LMERxLong<-ROIS$LMERxLong-ROIS$age ; ggplot(ROIS, aes(x=MAE_LMERxLong))+geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.6, fill="#0b6623")+geom_vline(xintercept = 0, size=2) + xlim(-10, 40) + theme_classic()

#########
### Prepare Slide 2 of Brain Age Gaps By Scan Sessions For All Processing Streams
#########

FSxCross<-ROIS[,c("sub","ses","age","MAE_FSxCross")] ; names(FSxCross)[4]<-"BrainAgeGap" ; FSxCross$STREAM<-"FSxCross" 
GAMxCross<-ROIS[,c("sub","ses","age","MAE_GAMxCross")] ; names(GAMxCross)[4]<-"BrainAgeGap" ; GAMxCross$STREAM<-"GAMxCross" 
LMERxCross<-ROIS[,c("sub","ses","age","MAE_LMERxCross")] ; names(LMERxCross)[4]<-"BrainAgeGap" ; LMERxCross$STREAM<-"LMERxCross" 
CROSS<-rbind(FSxCross,GAMxCross,LMERxCross) ; CROSS$ses<-as.factor(CROSS$ses) ; CROSS$STREAM<-as.factor(CROSS$STREAM)
FSxLong<-ROIS[,c("sub","ses","age","MAE_FSxLong")] ; names(FSxLong)[4]<-"BrainAgeGap" ; FSxLong$STREAM<-"FSxLong" 
GAMxLong<-ROIS[,c("sub","ses","age","MAE_GAMxLong")] ; names(GAMxLong)[4]<-"BrainAgeGap" ; GAMxLong$STREAM<-"GAMxLong" 
LMERxLong<-ROIS[,c("sub","ses","age","MAE_LMERxLong")] ; names(LMERxLong)[4]<-"BrainAgeGap" ; LMERxLong$STREAM<-"LMERxLong" 
Long<-rbind(FSxLong,GAMxLong,LMERxLong) ; Long$ses<-as.factor(Long$ses) ; Long$STREAM<-as.factor(Long$STREAM)
FINAL<-rbind(CROSS,Long)

ggplot(FINAL, aes(x=ses, y=BrainAgeGap,color=STREAM)) + geom_boxplot() + geom_abline(intercept=0,slope=0) + scale_color_manual(values=c("#de9f85","#F05E23","#bf3900","#81c795","#27a84c","#004713")) + theme_classic() + facet_wrap(~STREAM)

SUMMARY<-ROIS[,c("ses","MAE_FSxCross","MAE_LMERxCross")]
for (NUM in 1:5){
	print(paste0("CALCULATING SUMMARY STATS FOR SESSION:",NUM))
	ONE<-SUMMARY[which(SUMMARY$ses==NUM),] ; 
	print(mean(ONE$MAE_FSxCross)) ; print(sd(ONE$MAE_FSxCross)) 
	print(mean(ONE$MAE_LMERxCross)) ; print(sd(ONE$MAE_LMERxCross))
}	

TWO<-SUMMARY[which(SUMMARY$ses==2),] ; mean(TWO$MAE_FSxCross) ; sd(TWO$MAE_FSxCross) ; mean(TWO$MAE_LMERxCross) ; sd(TWO$MAE_LMERxCross) 
THREE<-SUMMARY[which(SUMMARY$ses==3),] ; mean(ONE$MAE_FSxCross) ; sd(ONE$MAE_FSxCross) ; mean(ONE$MAE_LMERxCross) ; sd(ONE$MAE_LMERxCross) 
FOUR<-SUMMARY[which(SUMMARY$ses==4),] ; mean(ONE$MAE_FSxCross) ; sd(ONE$MAE_FSxCross) ; mean(ONE$MAE_LMERxCross) ; sd(ONE$MAE_LMERxCross) 
FIVE<-SUMMARY[which(SUMMARY$ses==5),] ; mean(ONE$MAE_FSxCross) ; sd(ONE$MAE_FSxCross) ; mean(ONE$MAE_LMERxCross) ; sd(ONE$MAE_LMERxCross) 

#########
### Prepare Slide 3 of New Brain Age Figures with the best Harmonization Processing Stream
#########

REFINE<-FINAL[which(FINAL$STREAM == "FSxCross" | FINAL$STREAM == "LMERxCross"),]
FSxCross<-ROIS[,c("sub","ses","age","FSxCross","MAE_GAMxLong")] ; names(FSxCross)[4]<-"BrainAge" ; names(FSxCross)[5]<-"BrainAgeGap" ; FSxCross$STREAM<-"FSxCross" 
LMERxCross<-ROIS[,c("sub","ses","age","LMERxCross","MAE_GAMxLong")] ; names(LMERxCross)[4]<-"BrainAge" ; names(LMERxCross)[5]<-"BrainAgeGap" ; LMERxCross$STREAM<-"LMERxCross" 
HARM<-rbind(FSxCross,LMERxCross) ; HARM$ses<-as.numeric(HARM$ses) ; TRIO<-HARM[ which(HARM$ses < 4),] ; PRISMA<-HARM[ which(HARM$ses > 3),] ; HARM$ses<-as.factor(HARM$ses)

ggplot(HARM, aes(x=age,y=BrainAge,group=sub,color=STREAM)) + stat_smooth(method = "lm",se=F,size=0.6,alpha=0.5) + xlim(0,55) + ylim(0,55) + geom_abline(intercept=0,slope=1) + scale_color_manual(values=c("#de9f85","#bf3900")) + theme_classic() + facet_wrap(~STREAM)

ggplot(TRIO, aes(x=age,y=BrainAge,color=STREAM)) + stat_smooth(method = "lm",se=T,size=2,alpha=0.5,fullrange=T) + xlim(0,55) + ylim(0,55) + geom_abline(intercept=0,slope=1) + geom_point() + scale_color_manual(values=c("#de9f85","#bf3900")) + theme_classic() + facet_wrap(~STREAM)

ggplot(PRISMA, aes(x=age,y=BrainAge,color=STREAM)) + stat_smooth(method = "lm",se=T,size=2,alpha=0.5,fullrange=T) + xlim(0,55) + ylim(0,55) + geom_abline(intercept=0,slope=1) + geom_point() + scale_color_manual(values=c("#de9f85","#bf3900")) + theme_classic() + facet_wrap(~STREAM)
data$ManualQA<-as.factor(data$ManualQA)


geom_smooth(method = "lm", se = FALSE)

#########
### Prepare Slide 4 Showing Manual QA Exclusions Based On Raw T1w Scans
#########

data$BAG<-data$age
mu <- ddply(data, "SITE", summarise, grp.mean=mean(BAG))
ggplot(data, aes(x=BAG, color=SITE, fill=SITE)) +
	geom_histogram(aes(y=..density..), position="identity", alpha=0.80)+
	geom_density(alpha=0.2)+
	geom_vline(data=mu, aes(xintercept=grp.mean, color=SITE),linetype="longdash",size=1.9)+
	scale_color_manual(values=c("#000000","#bf3900", "#000000"))+
	scale_fill_manual(values=c("#000000","#bf3900", "#000000"))+
	theme_classic()

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
