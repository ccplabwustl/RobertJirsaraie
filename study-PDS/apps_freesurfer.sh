#!/usr/bin/env bash
###################

STUDY_DIR=/scratch/rjirsara/study-PDS
DIR_TOOLBOX=/scratch/rjirsara/RobertJirsaraie/toolbox
mkdir -p ${STUDY_DIR}/apps/freesurfer

module purge ; module load "freesurfer-5.3.0"

######
### Build Freesurfer License If Missing
######

FREESURFER_LICENSE=`echo $DIR_TOOLBOX/bids_apps/freesurfer/license_freesurfer.txt`
if [[ ! -f $FREESURFER_LICENSE && ! -z $DIR_TOOLBOX ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  "
	echo "`basename $FREESURFER_LICENSE` Not Found - Register For One Here: "
	echo "       https://surfer.nmr.mgh.harvard.edu/registration.html       "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  "
	printf "rjirsara@uci.edu
		40379 
		*CBT0GfF/00EU 
		FSP82KoWQu0tA \n" | grep -v local > $FREESURFER_LICENSE
fi

######
### Submit Jobs For Cross-Sectional Processing
######

SCRIPT_CROSS=${DIR_TOOLBOX}/bids_apps/freesurfer/pipeline_anatomical-cross.sh
if [[ -f $SCRIPT_CROSS && -d $STUDY_DIR ]] ; then
	for INPUT in `find ${STUDY_DIR}/bids -iname *_T1w.nii.gz | grep -v '_run-'` ; do
		SUB_IDS=`basename $INPUT | sed s@'_T1w.nii.gz'@''@g`
		SUB=`echo $SUB_IDS | cut -d '_' -f1 | cut -d '-' -f2`
		SES=`echo $SUB_IDS | cut -d '_' -f2 | cut -d '-' -f2`
		JOBNAME=`echo FSCROSS${SUB}x${SES}` ; JOBSTATUS=`qstat -u $USER | grep "${JOBNAME}\b" | awk {'print $10'}`
		if [[ ! -z "$JOBSTATUS" ]] ; then 
			echo ""
			echo "##########################################################"
			echo "#${SUB}x${SES} Is Currently Being Processed: ${JOBSTATUS} "
			echo "##########################################################"
		elif [[ -d  `echo $STUDY_DIR/apps/freesurfer/sub-${SUB}_ses-${SES}` ]] ; then
			echo ""
			echo "##########################################"
			echo "#${SUB}x${SES} Was Processed Successfully "
			echo "##########################################"
		else
			echo ""
			echo "####################################################################"
			echo "Submitting Freesurfer Cross-sectional Job For Subect: ${SUB}x${SES} "
			echo "####################################################################"
			cat $SCRIPT_CROSS \
				| sed s@'$1'@${STUDY_DIR}@g \
				| sed s@'$2'@${SUB}@g \
				| sed s@'$3'@${SES}@g > ${SCRIPT_CROSS}_GO
			qsub -N $JOBNAME ${SCRIPT_CROSS}_GO ; rm ${SCRIPT_CROSS}_GO
		fi
	done
fi

######
### Submit Jobs To Create Subject-Specific Base Templates
######

SCRIPT_BASE=${DIR_TOOLBOX}/bids_apps/freesurfer/pipeline_anatomical-base.sh
if [[ -f $SCRIPT_BASE && -d $STUDY_DIR ]] ; then
	for SUB in `find ${STUDY_DIR}/bids -maxdepth 1 | grep sub | sed s@'sub-'@'#'@g | cut -d '#' -f2` ; do
		NBID=`find ${STUDY_DIR}/bids/sub-${SUB} -iname *_T1w.nii.gz | grep -v _run- | wc -l`
		NFREE=`echo $STUDY_DIR/apps/freesurfer/sub-${SUB}_ses-* | tr ' ' '\n' | grep -v ERROR | grep -v long | wc -l`
		JOBNAME=`echo FSBASE${SUB}x${SES}` ; JOBSTATUS=`qstat -u $USER | grep "${JOBNAME}\b" | awk {'print $10'}`
		STATS_FILE=`echo $STUDY_DIR/apps/freesurfer/sub-${SUB}_base/stats/aseg.stats`
		if [[ ! -z "$JOBSTATUS" ]] ; then 
			echo ""
			echo "########################################################"
			echo "#${SUB}_base Is Currently Being Processed: ${JOBSTATUS} "
			echo "########################################################"
		elif [[ -f  $STATS_FILE ]] ; then
			echo ""
			echo "########################################"
			echo "#${SUB}_base Was Processed Successfully "
			echo "########################################"
		elif [[ $NBID != $NFREE ]] ; then 
			echo ""
			echo "##################################################################"
			echo "#${SUB}_base Does Not Have All Sessions Processed $NBID != $NFREE "
			echo "##################################################################"
		else
			echo ""
			echo "##################################################"
			echo "Submitting Freesurfer Base Job For Subect: ${SUB} "
			echo "##################################################"
			cat $SCRIPT_BASE \
				| sed s@'$1'@${STUDY_DIR}@g \
				| sed s@'$2'@${SUB}@g > ${SCRIPT_BASE}_GO
			qsub -N $JOBNAME ${SCRIPT_BASE}_GO ; rm ${SCRIPT_BASE}_GO
		fi
	done
fi

######
### Submit Jobs For Longitudinal Processing
######

SCRIPT_LONG=${DIR_TOOLBOX}/bids_apps/freesurfer/pipeline_anatomical-long.sh
if [[ -f $SCRIPT_LONG && -d $STUDY_DIR ]] ; then
	for SUB in `find ${STUDY_DIR}/apps/freesurfer -maxdepth 1 -printf "%f\n" | grep sub | grep -v _base | grep -v _long | head -n1` ; do
		SUBID=`echo $SUB | cut -d '_' -f1 | cut -d '-' -f2`
		SESID=`echo $SUB | cut -d '_' -f2 | cut -d '-' -f2`
		LONG_FILE=`echo ${STUDY_DIR}/apps/freesurfer/${SUB}_long/scripts/recon-all.log`
		BASE_FILE=`echo ${STUDY_DIR}/apps/freesurfer/sub-${SUBID}_base/scripts/recon-all.log`
		JOBNAME=`echo FSLONG${SUBID}x${SESID}` ; JOBSTATUS=`qstat -u $USER | grep "${JOBNAME}\b" | awk {'print $10'}`
		if [[ ! -z "$JOBSTATUS" ]] ; then 
			echo ""
			echo "########################################################"
			echo "#${SUB}_long Is Currently Being Processed: ${JOBSTATUS} "
			echo "########################################################"
		elif [[ -f ${LONG_FILE} ]] ; then
			echo ""
			echo "########################################"
			echo "#${SUB}_long Was Processed Successfully "
			echo "########################################"
		elif [[ ! -f "$BASE_FILE" ]] ; then 
			echo ""
			echo "#############################################################"
			echo "#${SUB} Does Not Have a Base Template Yet So Must Be Skipped "
			echo "#############################################################"
		else
			echo ""
			echo "##########################################################"
			echo "Submitting Freesurfer Longitudinal Job For Subect: ${SUB} "
			echo "##########################################################"
			cat $SCRIPT_LONG \
				| sed s@'$1'@${STUDY_DIR}@g \
				| sed s@'$2'@${SUBID}@g \
				| sed s@'$3'@${SESID}@g > ${SCRIPT_LONG}_GO
			qsub -N $JOBNAME ${SCRIPT_LONG}_GO ; rm ${SCRIPT_LONG}_GO
		fi
	done
fi

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
