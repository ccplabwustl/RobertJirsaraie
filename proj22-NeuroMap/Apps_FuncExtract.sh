#!/usr/bin/env bash
source ~/Anaconda_Env.sh
#####################

SUB=SUBID
ROOT=/scratch/jirsaraie
HCD=/ceph/intradb/archive/CinaB/CCF_HCD_STG
XCP=$ROOT/RobertJirsaraie/toolbox/bids_apps/dependencies/utils_xcp
SIMG=$ROOT/RobertJirsaraie/toolbox/bids_apps/dependencies/xcpengine.simg
ATLAS=$ROOT/RobertJirsaraie/toolbox/bids_apps/atlases/atl-Gordon_T1w_Parcels.nii.gz
LABELS=$ROOT/RobertJirsaraie/toolbox/bids_apps/atlases/atl-Gordon_Labels_Parcels.1D
PREPROC=$ROOT/study-HCD/bids/$SUB/func/${SUB}_rfMRI_REST_hp0_clean_mgtr_hpss.nii.gz
CENSORED=$ROOT/study-HCD/bids/$SUB/func/${SUB}_fully_denoised.nii.gz

#####
### Aggregate the XCP Utilities/Atlas
#####

if [[ ! -d $XCP || ! -f $SIMG ]] ; then 
	singularity build ${SIMG} docker://pennbbl/xcpengine:latest
	mkdir -p $XCP ; mv xcpEngine/utils/* $XCP ; rm -rf xcpEngine
fi

#####
### Temporally Censor the High Motion TRs
#####

if [[ -f $PREPROC && ! -f $CENSORED ]] ; then
	echo "${SUB} - CENSORED fMRI" ; cd `dirname $CENSORED`
	fslsplit $PREPROC -t ; rm -rf FDfilt4_censored.txt
	COUNT=`fslinfo $PREPROC | grep dim4 | head -n1 | awk {'print $2'}`
	for INDEX in `seq 1 1 $COUNT` ; do 
		FD=`cat FDfilt4.txt | head -n${INDEX} | tail -n1`
		VOLUME=`ls vol*.nii.gz | head -n${INDEX} | tail -n1`
		if [[ $FD < 0.2 ]] ; then
			echo $FD >> FDfilt4_censored.txt
		else
			rm $VOLUME
		fi
	done
	fslmerge -t $CENSORED vol*.nii.gz ; rm vol*.nii.gz $PREPROC
fi

#####
### Compute Functional Connectivity Matrix
#####

if [[ ! -f $ROOT/study-HCD/apps/xcp-fcon/${SUB}/fcon_parcels.tsv ]] ; then
	echo "${SUB} - FCON Module" ; mkdir -p $ROOT/study-HCD/apps/xcp-fcon/${SUB}
	TIME=$ROOT/study-HCD/apps/xcp-fcon/${SUB}/timeseries.1D
	singularity exec -B ${ROOT} ${SIMG} \
		/xcpEngine/utils/roi2ts.R \
		-i $CENSORED \
		-r $ATLAS > $TIME
	NET=$ROOT/study-HCD/apps/xcp-fcon/${SUB}/network.txt
	singularity exec -B ${ROOT} ${SIMG} \
		/xcpEngine/utils/ts2adjmat.R \
		-t $TIME > $NET
	MISSING=$ROOT/study-HCD/apps/xcp-fcon/${SUB}/problematic.1D
	singularity exec -B ${ROOT} ${SIMG} \
		/xcpEngine/utils/missingIdx.R \
		-i $NET > $MISSING
	find $ROOT/study-HCD/apps/xcp-fcon/${SUB} -type f -empty -print -delete
	FCON_PARCELS=$ROOT/study-HCD/apps/xcp-fcon/${SUB}/fcon_parcels.tsv 
	singularity exec -B ${ROOT} ${SIMG} \
		/xcpEngine/utils/adjmat2pajek.R \
		-a $NET \
		-t N | tr ' ' '\t' | tail -n+3 > $FCON_PARCELS ; rm $NET
fi

#####
### Compute Amplitude of Low Frequency Fluctuations Maps
#####

if [[ ! -f $ROOT/study-HCD/apps/xcp-alff/${SUB}/alff_quantifyAtlas.csv ]] ; then 
	echo "${SUB} - ALFF Module" ; mkdir -p $ROOT/study-HCD/apps/xcp-alff/${SUB}
	fslpspec $CENSORED ${SUB}_fslpspec.nii.gz
	fslmaths  ${SUB}_fslpspec.nii.gz -sqrt ${SUB}_fslpspec_sqrt.nii.gz
	fslroi ${SUB}_fslpspec_sqrt.nii.gz ${SUB}_fslpspec_sqrt_lowfreq.nii.gz 2 22
	ALFF=$ROOT/study-HCD/apps/xcp-alff/${SUB}/alff.nii.gz
	fslmaths ${SUB}_fslpspec_sqrt_lowfreq.nii.gz -Tmean -mul 22 $ALFF
	ZALFF=$ROOT/study-HCD/apps/xcp-alff/${SUB}/Zalff.nii.gz
	SD=`fslstats $ALFF -S` ; MEAN=`fslstats $ALFF -M` 
	fslmaths $ALFF -sub $MEAN -div $SD $ZALFF
	rm -rf ${SUB}_fslpspec*.nii.gz
	for SCALE in "alff" "Zalff" ; do
		QUANT=$ROOT/study-HCD/apps/xcp-alff/${SUB}/${SCALE}_quantifyAtlas.csv 
		singularity exec -B ${ROOT} ${SIMG} \
			/xcpEngine/utils/quantifyAtlas \
			-v $ROOT/study-HCD/apps/xcp-alff/${SUB}/${SCALE}.nii.gz \
			-a $ATLAS \
			-n $SCALE \
			-p $SUB \
			-s mean \
			-r $LABELS \
			-o $QUANT
		rm $ROOT/study-HCD/apps/xcp-alff/${SUB}/${SCALE}*.1D
	done
fi

#####
### Compute Extract Regional Homogeneity Maps
#####

if [[ ! -f $ROOT/study-HCD/apps/xcp-reho/${SUB}/reho_quantifyAtlas.csv ]] ; then 
	echo "${SUB} - Reho Module" ; mkdir -p $ROOT/study-HCD/apps/xcp-reho/${SUB}
	REHO=$ROOT/study-HCD/apps/xcp-reho/${SUB}/reho.nii.gz
	3dReHo -prefix $REHO -inset $CENSORED -nneigh 27
	ZREHO=$ROOT/study-HCD/apps/xcp-reho/${SUB}/Zreho.nii.gz
	MEAN=`fslstats $REHO -M` ; SD=`fslstats $REHO -S`
	fslmaths $REHO -sub $MEAN -div $SD $ZREHO
	for SCALE in "reho" "Zreho" ; do
		QUANT=$ROOT/study-HCD/apps/xcp-reho/${SUB}/${SCALE}_quantifyAtlas.csv 
		singularity exec -B ${ROOT} ${SIMG} \
			/xcpEngine/utils/quantifyAtlas \
			-v $ROOT/study-HCD/apps/xcp-reho/${SUB}/${SCALE}.nii.gz \
			-a $ATLAS \
			-n $SCALE \
			-p $SUB \
			-s mean \
			-r $LABELS \
			-o $QUANT
		rm $ROOT/study-HCD/apps/xcp-reho/${SUB}/${SCALE}*.1D
	done
fi

#####
### Compute Network Integrity Maps
#####

if [[ ! -f $ROOT/study-HCD/apps/xcp-dualreg/${SUB}/stats/node-1_cohesion.txt ]] ; then
	echo "${SUB} - DualReg Module"
	DRPath=$ROOT/study-HCD/apps/xcp-dualreg 
	ATL4D=$DRPath/atlas/atl-Gordon_T1w_Parcels4D.nii.gz
	NUM_NODES=`fslstats $ATLAS -R | cut -d ' ' -f2 | cut -d '.' -f1`
	if [[ ! -f $ATL4D ]] ; then
		mkdir -p $DRPath/atlas
		for NODE in `seq 1 1 $NUM_NODES` ; do
			fslmaths $ATLAS -thr $NODE -uthr $NODE -div $NODE $DRPath/atlas/atl_node-${NODE}.nii.gz
		done
		fslmerge -t $ATL4D $DRPath/atlas/atl_node-*.nii.gz
	fi
	rm -rf $DRPath/${SUB} ; mkdir -p $DRPath/${SUB}/stats
	dual_regression $ATL4D 1 -1 5000 $ROOT/study-HCD/apps/xcp-dualreg/${SUB} $CENSORED
	for NODE in `seq 1 1 $NUM_NODES` ; do
		BETA_MAP=`ls $DRPath/${SUB}/dr_stage2_ic*.nii.gz | head -n${NODE} | tail -n1`
		BETA_MAP_THR=`echo $BETA_MAP | sed s@'stage2'@'stage2-thr'@g`
		fslmaths $BETA_MAP -mul $DRPath/atlas/atl_node-${NODE}.nii.gz $BETA_MAP_THR
		fslstats $BETA_MAP_THR -M > $DRPath/${SUB}/stats/node-${NODE}_cohesion.txt
		fslstats $BETA_MAP -M > $DRPath/${SUB}/stats/node-${NODE}_integration.txt	
	done
fi

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######