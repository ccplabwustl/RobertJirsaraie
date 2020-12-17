#!/usr/bin/env Rscript
######################

print("Reading Arguments")

covaPath <- "/dfs2/yassalab/rjirsara/ConteCenterScripts/Conte-One/analyses/IntraFlux/n138_IntraFlux.mods/NetStr.csv" 
inputPath <- "/dfs2/yassalab/rjirsara/ConteCenterScripts/Conte-One/analyses/IntraFlux/n138_IntraFlux.mods/LONG.csv"
OutDirRoot <- "/dfs2/yassalab/rjirsara/ConteCenterScripts/Conte-One/analyses/IntraFlux/n138_IntraFlux.mods"
DATASETS<-list.files("/dfs2/yassalab/rjirsara/ConteCenterScripts/Conte-One/analyses/IntraFlux",pattern="Aggregate",full.names=TRUE)
DATASETS<-DATASETS[grepl("Longitudinal",DATASETS)]
m1 <- "~AgeAtScan+Gender+FD_MEAN+PreMood_Lvl"
m2 <- "~AgeAtScan+Gender+FD_MEAN+PreMood_Ent"
m3 <- "~AgeAtScan+Gender+FD_MEAN+PreMood_Lvl+PreMood_Ent"
m4 <- "~AgeAtScan+Gender+FD_MEAN+PreMood_Lvl+Gender*PreMood_Ent"
m5 <- "~AgeAtScan+Gender+FD_MEAN+PreMood_Lvl+AgeAtScan*PreMood_Ent"
m6 <- "~AgeAtScan+Gender+FD_MEAN+PreMood_Lvl+TASK*PreMood_Ent"
MODELS <- list(m1,m2,m3,m4,m5,m6)
randomFormula <- "~(1|sub)"

print("Reading Custom Outputting Options")
SubOutDir=FALSE
SortOutFile=FALSE
SaveDatasets=FALSE

#####################
### Load Packages ###
#####################

print("Loading Libraries")

suppressMessages(require(ggplot2))
suppressMessages(require(base))
suppressMessages(require(reshape2))
suppressMessages(require(nlme))
suppressMessages(require(lme4))
suppressMessages(require(gamm4))
suppressMessages(require(stats))
suppressMessages(require(knitr))
suppressMessages(require(mgcv))
suppressMessages(require(plyr))
suppressMessages(require(oro.nifti))
suppressMessages(require(parallel))
suppressMessages(require(optparse))
suppressMessages(require(fslr))
suppressMessages(require(voxel))
suppressMessages(require(stringr))

##############################################
### Load and Prepare Datasets for Analyses ###
##############################################

for (covsFormula in MODELS){

	print(covsFormula)
	print("Loading Covariates Dataset")
	covaData<-read.csv(DATASETS)
	covaData<-covaData[,c("sub","AgeAtScan","Gender","FD_MEAN","scl.CDI_MD","PreMood_Lvl","PreMood_Ent","TASK")]
	if (any(colnames(covaData) == "X")) {
		covaData$X<-NULL
	}

	print("Loading Responce Dataset")
	inputData <- read.csv(DATASETS)
	inputData<-inputData[grepl("_MEAN",names(inputData))]
	inputData<-inputData[grepl("COMP",names(inputData))]
	if (any(colnames(inputData) == "X")) {
		inputData$X<-NULL
	}

	print("Remove Empty Columns with No Variance")
	completecols<-which(colSums(inputData,na.rm=T)!=0)
	inputData<-inputData[,c(completecols)]

	print("Merging Datasets")
	dataSubj <- cbind(covaData, inputData)
	dataSubj<-dataSubj[-c(which(dataSubj$FD_MEAN > 0.5)),]

	print("Remove Missing Observations")
	Predictors<-gsub("~", "",covsFormula)
	Predictors<-gsub("\\*", "+",Predictors)
	Predictors<-strsplit(Predictors, "+", fixed = TRUE)
	MaxPredictors<-dim(as.data.frame(Predictors))[1]
	covaData$Gender<-as.factor(covaData$Gender)

	print("Ensure Class Type Of Predictors")
	for (x in 1:MaxPredictors){
		var<-as.character(as.data.frame(Predictors)[x,1])
		if (startsWith(var, 's(')){ 
			var<-strsplit(var, ",", fixed = TRUE)
			var<-as.character(as.data.frame(var)[1,1])
			var<-substring(var, 3)
		}
		if (var %in% colnames(covaData)){
			covaData<-covaData[complete.cases(covaData[,var]),]
			classtype<-class(covaData[,var])
			print(paste("########################################"))
			print(paste("⚡⚡⚡",var,"class type is", classtype,"⚡⚡⚡"))
			print(paste("########################################"))
		} else {
			print(paste("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
			print(paste("⚡⚡⚡",var,"NOT FOUND - EXITING SCRIPT ⚡⚡⚡"))
			print(paste("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"))
			exit()
		}
	}

#################################################
### Define Cores & Multiple Correction Method ###
#################################################

	ncores <- 5
	pAdjustMethod <- "BH"
	methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
	if (!any(pAdjustMethod == methods)) {
		print("p.adjust.method is not a valid one, reverting to 'none'")
		pAdjustMethod <- "none"
	}

########################################################
### Defining and Executing GAMM4 Mixed Effect Models ###
########################################################

	print("Analyzing Dataset")
	model.formula <- mclapply((dim(covaData)[2] + 1):dim(dataSubj)[2], function(x) { 
		as.formula(paste(paste0("dataSubj[,",x,"]"), covsFormula, sep="")) 
	}, mc.cores=ncores)

	print("Executing Models")
	m <- mclapply(model.formula, function(x) {
		ANALYZE <- gamm4(formula = x, random=as.formula(randomFormula), data=dataSubj, REML=T)$gam
		summary <- summary(ANALYZE)
		residuals <- ANALYZE$residuals
		missing <- as.numeric(ANALYZE$na.action)
		return(list(summary,residuals, missing))
	}, mc.cores=ncores)

#################################
### Creating Output Directory ###
#################################

	print("Creating Analysis Directory")

	subjDataOut <- strsplit(covaPath, ".csv")[[1]][[1]]
	subjDataOut <- strsplit(subjDataOut, "/")[[1]][[length(subjDataOut <- strsplit(subjDataOut, "/")[[1]])]]
	inputPathOut <- strsplit(inputPath, ".csv")[[1]][[1]]
	inputPathOut <- strsplit(inputPathOut, "/")[[1]][[length(inputPathOut <- strsplit(inputPathOut, "/")[[1]])]]
	OutDir <-  paste0(OutDirRoot,"/COVA-",subjDataOut,"_RESP-",inputPathOut)
	suppressWarnings(dir.create(OutDir, recursive = TRUE))
	setwd(OutDir)

	if (SubOutDir != 'FALSE') {
		print("Adding Custom Sub-Directory")
		OutDir<-paste0(OutDir,"/",SubOutDir)
		suppressWarnings(dir.create(OutDir))
		setwd(OutDir)
	}

	print("Creating Output File Name")

	outName <- gsub("~", "", covsFormula)
	outName <- gsub(" ", "", outName)
	outName <- gsub("\\+","+",outName)
	outName <- gsub("\\(","",outName)
	outName <- gsub("\\)","",outName)
	outName <- gsub(",","",outName)
	outName <- gsub("\\.","",outName)
	outName <- gsub("=","",outName)
	outName <- gsub("\\*","and",outName)
	outName <- gsub(":","and",outName)
	random <- gsub("~", "", randomFormula)
	random <- gsub("\\(", "", random)
	random <- gsub("\\)", "", random)
	random <- gsub("\\|", "", random)
	if (grepl("s*,*k=", covsFormula)){
		outsubDir <- paste0("n",dim(dataSubj)[1],"_gamm4_",outName,"_random_",random)
		outsubDir<-paste(OutDir,outsubDir,sep="/")
	} else {
		outsubDir <- paste0("n",dim(dataSubj)[1],"_lme4_",outName,"_random_",random)
		outsubDir<-paste(OutDir,outsubDir,sep="/")
	}

############################################
### Save Processed Datasets if Specified ###
############################################

	if (SaveDatasets == 'TRUE') {
		print("Saving Processed Datasets")
		DataOutDir<-str_replace_all(OutDir, "/Results/", "/Data/")
		suppressWarnings(dir.create(DataOutDir, recursive = TRUE))
		FileName<-paste(basename(outsubDir),"Predictors.csv",sep="_")
		write.csv(covaData,paste(DataOutDir,FileName,sep="/"))
		FileName<-paste(basename(outsubDir),"Responces.csv",sep="_")
		write.csv(inputData,paste(DataOutDir,FileName,sep="/"))
		FileName<-paste(basename(outsubDir),"Merged.csv",sep="_")
		write.csv(dataSubj,paste(DataOutDir,FileName,sep="/"))
		Sys.chmod(list.files(DataOutDir, full.names=TRUE), "775", use_umask = FALSE)
	}

##########################################
### Generating Table of Summary Output ###
##########################################

	print("Generating parameters")

	m <- mclapply(m, function(x) {
		x[[1]]
	}, mc.cores=ncores)
	
	length.names.p <- length(rownames(m[[1]]$p.table))
	output <- as.data.frame(matrix(NA, 
			nrow = length((dim(covaData)[2] + 1):dim(dataSubj)[2]), 
			ncol= (1+3*length.names.p)))
	names(output)[1] <- "names"

	for (i in 1:length.names.p){
		dep.val <- rownames(m[[1]]$p.table)[i]
		names(output)[2 + (i-1)*3 ] <- paste0("tval.",dep.val)
		names(output)[3 + (i-1)*3 ] <- paste0("pval.",dep.val)
		names(output)[4 + (i-1)*3 ] <- paste0("pval.",pAdjustMethod,dep.val)
		val.tp <- t(mcmapply(function(x) {
			x$p.table[which(rownames(x$p.table) == dep.val), 3:4]
		}, m, mc.cores=ncores))
		output[,(2 + (i-1)*3):(3 + (i-1)*3)] <- val.tp
		output[,(4 + (i-1)*3)] <- p.adjust(output[,(3 + (i-1)*3)], pAdjustMethod)
	}

	output$names <- names(dataSubj)[(dim(covaData)[2] + 1):dim(dataSubj)[2]]
	p.output <- output

####################################
### Save Final Output of Results ###
####################################

	print("Saving Final Results File")

	if (is.null(m[[1]]$s.table)){  
		if (SortOutFile != 'FALSE'){
			COLNAME<-grep(SortOutFile, names(p.output), value = TRUE)[1]
			COLNUM<-which(names(p.output) == COLNAME )
			p.output<-p.output[order(p.output[,COLNUM] ),]
		}
		outsubDir<-str_replace_all(outsubDir, "_gamm4_", "_lme4_")
		write.csv(p.output, paste0(outsubDir,".csv"), row.names=F)
		Sys.chmod(list.files(dirname(outsubDir), full.names=TRUE), "775", use_umask = FALSE)
	} else {
		length.names.s <- length(rownames(m[[1]]$s.table))
		output <- as.data.frame(matrix(NA, 
				nrow = length((dim(covaData)[2] + 1):dim(dataSubj)[2]), 
				ncol= (1+2*length.names.s)))
		names(output)[1] <- "names"
		for (i in 1:length.names.s){
			dep.val <- rownames(m[[1]]$s.table)[i]
			names(output)[2 + (i-1)*2 ] <- paste0("pval.",dep.val)
			names(output)[3 + (i-1)*2 ] <- paste0("pval.",pAdjustMethod,dep.val)
			val.tp <- mcmapply(function(x) {
				x$s.table[which(rownames(x$s.table) == dep.val), 4]
			}, m, mc.cores=ncores)
   			output[,(2 + (i-1)*2)] <- val.tp
			output[,(3 + (i-1)*2)] <- p.adjust(output[,(2 + (i-1)*2)], pAdjustMethod) 
			output$names <- names(dataSubj)[(dim(covaData)[2] + 1):dim(dataSubj)[2]]
			output <- merge(p.output, output, by="names")
			if (SortOutFile != 'FALSE') {
				COLNAME<-grep(SortOutFile, names(p.output), value = TRUE)[1]
				COLNUM<-which(names(p.output) == COLNAME )
				p.output<-p.output[order(p.output[,COLNUM] ),]
			}
			write.csv(output, paste0(outsubDir,".csv"), row.names=F)
			Sys.chmod(list.files(dirname(outsubDir), full.names=TRUE), "775", use_umask = FALSE)
		}
	}

	print("Script Ran Succesfully")

}

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
