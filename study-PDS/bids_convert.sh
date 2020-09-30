#!/bin/bash
###########

DIR_PROJECT=/scratch/rjirsara/study-PDS
DIR_GITHUB=/scratch/rjirsara/RobertJirsaraie
mkdir -p $DIR_PROJECT/dicoms $DIR_PROJECT/bids $DIR_PROJECT/audits

##############################################
### Convert dcm2bids for the first 3 Waves ###
##############################################

SCRIPT_DCM2BIDS=${DIR_GITHUB}/toolbox/dicom_management/bids_convertDCM.sh
FILE_CONFIG=${DIR_GITHUB}/study-PDS/bids_config_waves123.json
OPT_LONGITUDINAL=sub_ses
OPT_RM_DICOMS=TRUE

if [[ -f $SCRIPT_DCM2BIDS && -f $FILE_CONFIG ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡"
	echo "Coverting Data into BIDs Format - dcm2bids "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡"
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
	${SCRIPT_DCM2BIDS} ${FILE_CONFIG} ${DIR_PROJECT}/dicoms ${OPT_LONGITUDINAL} ${OPT_RM_DICOMS}
fi

########################################################
### Modify MetaData to Adhere to BIDS Specifications ###
########################################################

SCRIPT_EDIT_META=${DIR_GITHUB}/dicom_management/bids_conversion/MetaData_Modification.py
OPT_GEN_FMAP_FUNC=FALSE
OPT_GEN_FMAP_DWI=FALSE
OPT_ADD_DEFAULT_ST=TRUE

if [[ -f $SCRIPT_EDIT_META ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	echo "Editing Meta Data to Adhere to BIDS Specifications "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	python $SCRIPT_EDIT_META $DIR_PROJECT $OPT_GEN_FMAP_FUNC $OPT_GEN_FMAP_DWI $OPT_ADD_DEFAULT_ST
fi

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

#################################
### Upload Events To Flywheel ###
#################################

SCRIPT_FW_UPLOAD=${DIR_GITHUB}/dicom_management/flywheel/Flywheel_Upload.sh
DIR_FLYWHEEL=yassalab/Conte-Two-UCI
DIRNAMExTYPE="DICOMSxdicom PARRECxparrec NIFTIxfolder"

if [[ -f $SCRIPT_FW_UPLOAD ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	echo "Uploading Events To Subject-Level Directories In Flyweel "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	$SCRIPT_FW_UPLOAD $DIR_PROJECT $DIR_FLYWHEEL $DIRNAMExTYPE $FLYWHEEL_API_TOKEN
fi

########################################
### Uploading Processed Data To Zion ###
########################################

SCRIPT_ZION_UPLOAD=${DIR_GITHUB}/dicom_management/zion/Zion_Upload.exp
DIR_HPC=${DIR_PROJECT}/bids
DIR_ZION=/tmp/yassamri/Conte_Center/Processed_Data/rjirara
ZION_USERNAME=`whoami`

if [[ -f $SCRIPT_ZION_UPLOAD && ! -z $DIR_ZION && ! -z $ZION_PASSWORD ]] ; then
	echo ""
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	echo "Uploading Processed Data To Zion "
	echo "⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  #  ⚡  "
	$SCRIPT_ZION_UPLOAD $DIR_HPC $DIR_ZION $ZION_USERNAME $ZION_PASSWORD
fi

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
######         ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡       ######
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
