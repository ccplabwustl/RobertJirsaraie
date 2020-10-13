#!/bin/bash
###########

DIR_PROJECT=/scratch/rjirsara/study-PDS
DIR_GITHUB=/scratch/rjirsara/RobertJirsaraie
mkdir -p $DIR_PROJECT/dicoms $DIR_PROJECT/bids $DIR_PROJECT/audits

###########################################################
### Download Dicoms form the cnda for the first 3 Waves ###
###########################################################

SCRIPT_CNDA=${DIR_GITHUB}/toolbox/dicom_management/cnda_download.sh

if [[ -f $SCRIPT_CNDA ]] ; then
	for ZIP in `echo 105046` ; do
		ZIP_DOWNLOAD=`echo ${USER}-$(date "+%Y%m%d")_${ZIP}`
		echo ""
		echo "############################################################"
		echo "Submitting Download Job ${ZIP_DOWNLOAD} To Aggregate Dicoms "
		echo "############################################################"
	 	cat $SCRIPT_CNDA | \
			sed s@'$1'@"${DIR_PROJECT}/dicoms"@g | \
			sed s@'$2'@"${ZIP_DOWNLOAD}"@g > ${ZIP_DOWNLOAD}.sh
		qsub -N 'CNDA' ${ZIP_DOWNLOAD}.sh ; rm ${ZIP_DOWNLOAD}.sh
	done
fi

##############################################
### Convert dcm2bids for the first 3 Waves ###
##############################################

SCRIPT_DCM2BIDS=${DIR_GITHUB}/toolbox/dicom_management/bids_convertDCM.sh
FILE_CONFIG=${DIR_GITHUB}/study-PDS/bids_config123.json
OPT_LONGITUDINAL=sub_ses
OPT_RM_DICOMS=TRUE

if [[ -f $SCRIPT_DCM2BIDS && -f $FILE_CONFIG ]] ; then
	if [[ -d `echo $DIR_PROJECT/dicoms/* | tr ' ' '\n' | grep _L | head -n1` ]] ; then
		SUBJECTS=`echo $DIR_PROJECT/dicoms/* | sed s@"$DIR_PROJECT/dicoms/"@''@g`
		for SUB in `echo $SUBJECTS | tr ' ' '\n' | cut -d '_' -f2 | sed s@'L'@@g | sort | uniq` ; do
			INDEX=0
			for SES in $DIR_PROJECT/dicoms/*L${SUB} ; do
				INDEX=$(($INDEX+1))
				mv $SES `echo ${SES}_${INDEX} | cut -d 'L' -f2`
				mkdir ${SES}_${INDEX}/dicoms
				mv ${SES}_${INDEX}/* ${SES}_${INDEX}/dicoms
			done
		done
	fi
	for DIR_SUBJECT in `echo ${DIR_PROJECT}/dicoms/*/dicoms | tr ' ' '\n' | sed s@"/dicoms$"@""@g` ; do
		JOBNAME=`echo BIDS_$(basename ${DIR_SUBJECT})| cut -c1-10`
		JOBSTATUS=`qstat -u $USER | grep "${JOBNAME}\b" | awk {'print $10'}`
		if [ -d $DIR_SUBJECT/sub-* ] ; then
			echo ""
			echo "#######################################"
			echo "#${JOBNAME} has already been processed "
			echo "#######################################"
		elif [[ ! -z "$JOBSTATUS" &&  $JOBSTATUS != 'C' ]] ; then 
			echo ""
			echo "#######################################################"
			echo "#${JOBNAME} Is Currently Being Processed: ${JOBSTATUS} "
			echo "#######################################################"
		else
			echo ""
			echo "##################################################"
			echo "Submitting DCM2BIDS Job ${JOBNAME} For Converting "
			echo "##################################################"
			cat ${SCRIPT_DCM2BIDS} | \
				sed s@'$1$'@"${DIR_SUBJECT}"@g | \
				sed s@'$2$'@"${FILE_CONFIG}"@g | \
				sed s@'$3$'@"${OPT_LONGITUDINAL}"@g | \
				sed s@'$4$'@"${OPT_RM_DICOMS}"@g > $DIR_SUBJECT/bids_convertDCM.sh
			qsub -N ${JOBNAME} $DIR_SUBJECT/bids_convertDCM.sh
		fi
	done	
fi

########################################################
### Modify MetaData to Adhere to BIDS Specifications ###
########################################################

SCRIPT_EDIT_META=${DIR_GITHUB}/toolbox/dicom_management/bids_modifyMETA.py
OPT_GEN_FMAP_FUNC=FALSE
OPT_GEN_FMAP_DWI=FALSE
OPT_ADD_DEFAULT_ST=FALSE

if [[ -f $SCRIPT_EDIT_META ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	echo "Editing Meta Data to Adhere to BIDS Specifications "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	python $SCRIPT_EDIT_META $DIR_PROJECT $OPT_GEN_FMAP_FUNC $OPT_GEN_FMAP_DWI $OPT_ADD_DEFAULT_ST
fi

<<SKIP
#######################################
### Generate BIDS Validation Report ###
#######################################

SCRIPT_GEN_REPORT=${DIR_GITHUB}/dicom_management/bids_conversion/Validation_Report.sh

if [[ -f $SCRIPT_GEN_REPORT ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  "
	echo "Generating BIDs Validation Report to Assess Standardization "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  "
	$SCRIPT_GEN_REPORT $DIR_PROJECT
fi

SKIP
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
