#!/usr/bin/env python
#PBS -o /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stdout
#PBS -e /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stderr
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=1:00:00
#PBS -N COMBAT_LME4
#PBS -q dque
#PBS -j oe 
#source ~/ConfigEnv.sh
#####################
#install.packages("devtools")
#devtools::install_github("jcbeer/longCombat")

library("longCombat")
PROJECT_DIR="/scratch/rjirsara/proj20-BrainAgeEval"

######
### Prepare Datasets to be Harmonization 
######

COVARS<-read.csv(paste0(PROJECT_DIR,"/analysis/n834_DataFreeze_20201213.csv"))
CROSS<-read.csv(paste0(PROJECT_DIR,"/analysis/n829_FS-Cross_20201231.csv"))
LONG<-read.csv(paste0(PROJECT_DIR,"/analysis/n828_FS-Long_20201231.csv"))
COVARS$SCANNER<-0
COVARS[which(COVARS$ses < 4),"SCANNER"]<-1
COVARS[which(COVARS$ses > 3),"SCANNER"]<-2
COVARS<-COVARS[,c("sub","ses","age","sex","SCANNER")]
COVARS$sex<-as.factor(COVARS$sex) ; COVARS$SCANNER<-as.factor(COVARS$SCANNER) 
LONG<-merge(COVARS,LONG,by=c("sub","ses","sex")) ; CROSS<-merge(COVARS,CROSS,by=c("sub","ses","sex"))

######
### Run Longitudinal Harmonization on Both Freesurfer Processing Streams
######

COMBAT_CROSS=longCombat(
	idvar="sub",
	timevar="ses",
	batchvar="SCANNER",
	features=names(CROSS)[6:ncol(CROSS)],
	formula="age+sex",
	ranef="(1|sub)",
	data=CROSS,
	niter = 30,
	method = "REML",
	verbose = TRUE
)

COMBAT_LONG=longCombat(
	idvar="sub",
	timevar="ses",
	batchvar="SCANNER",
	features=names(LONG)[6:ncol(LONG)],
	formula="age+sex",
	ranef="(1|sub)",
	data=LONG,
	niter = 30,
	method = "REML",
	verbose = TRUE
)

FINAL_LONG<-merge(LONG[,c("sub","ses","age","sex")],COMBAT_LONG$data_combat,by=c("sub","ses"))
FINAL_CROSS<-merge(CROSS[,c("sub","ses","age","sex")],COMBAT_CROSS$data_combat,by=c("sub","ses"))
write.csv(FINAL_CROSS,paste0(PROJECT_DIR,"/analysis/n829_COMBATxLMER-Cross_20210101.csv"),row.names=F)
write.csv(FINAL_LONG,paste0(PROJECT_DIR,"/analysis/n828_COMBATxLMER-Long_20210101.csv"),row.names=F)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
