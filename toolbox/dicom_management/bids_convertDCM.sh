#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=4:00:00
#PBS -q dque
#PBS -o ~
#PBS -e ~
source ~/ConfigEnv.sh
############

DIR_SUBJECT=$1
FILE_CONFIG=$2
OPT_LONGITUDINAL=$3
OPT_RM_DICOMS=$4

####################################
### Covert Dicoms To BIDs Format ###
####################################

for DICOMS in $DIR_SUBJECT ; do
	SUBID=`basename $DICOMS`
	DIR_PROJECT=`dirname $DICOMS`
	echo "Now Coverting Dicoms to BIDS for ${SUBID}"
	if [[ -z $OPT_LONGITUDINAL || $OPT_LONGITUDINAL == "FALSE" ]] ; then
		dcm2bids -d ${DIR_PROJECT}/${SUBID}/dicoms \
			-p ${SUBID} \
			-c ${FILE_CONFIG} \
			-o ${DIR_PROJECT}/${SUBID}/sub-${SUBID} \
			--forceDcm2niix \
			--clobber 
	elif [[ $OPT_LONGITUDINAL == *sub* && $OPT_LONGITUDINAL == *ses* ]] ; then
		delimiter=`echo $OPT_LONGITUDINAL | sed s@'^.*sub'@''@g |  sed s@"ses.*$"@@g`
		sub=`echo $SUBID | cut -d "${delimiter}" -f1`  
		ses=`echo $SUBID | cut -d "${delimiter}" -f2`
		dcm2bids -d ${DIR_PROJECT}/${SUBID}/dicoms \
			-p ${sub} \
			-s ${ses} \
			-c ${FILE_CONFIG} \
			-o ${DIR_PROJECT}/${SUBID}/sub-${SUBID} \
			--forceDcm2niix \
			--clobber 
	else
		dcm2bids -d ${DIR_PROJECT}/${SUBID}/dicoms \
			-p ${SUBID} \
			-s ${OPT_LONGITUDINAL} \
			-c ${FILE_CONFIG} \
			-o ${DIR_PROJECT}/${SUBID}/sub-${SUBID} \
			--forceDcm2niix \
			--clobber 
	fi

######################################################
### Move Field Maps With the Wrong Phase Direction ###
######################################################

	echo "Checking Phase Encoding Directions of Field Maps For ${SUBID}"
	DIR_LOCAL_FMAP=`find ${DIR_PROJECT}/${SUBID}/sub-${SUBID} -iname fmap`
	if [[ -d ${DIR_LOCAL_FMAP} && -f `echo $DIR_LOCAL_FMAP/*_dir-* | cut -d ' ' -f1` ]] ; then
		PHASE_DIRECTIONS="AP_j- PA_j"
		for DIRECTION in $PHASE_DIRECTIONS ; do
			PHASE=`echo ${DIRECTION} | cut -d '_' -f1`
			DIRECTION=`echo ${DIRECTION} | cut -d '_' -f2`
			for JSON in `echo $DIR_LOCAL_FMAP/*[da][ic][rq]-${PHASE}*.json` ; do
				if [[ -f $JSON ]] ; then
					FILE_PHASE=`cat $JSON | grep "PhaseEncodingDir" | grep -v "InPlane" | awk '{$1=$1;print}' | cut -d ' ' -f2 | sed s@'"'@@g | sed s@','@@g`
					if [[ $FILE_PHASE != ${DIRECTION} ]] ; then
						PROBLEMATIC_FILES=`echo $JSON | sed s@'json'@'*'@g`
						mkdir -p ${DIR_PROJECT}/${SUBID}/sub-${SUBID}/tmp_wrongphase
						mv $PROBLEMATIC_FILES ${DIR_PROJECT}/${SUBID}/sub-${SUBID}/tmp_wrongphase
					elif [[ $JSON == *"_run-"* ]] ; then
						FILEINDEX=`basename ${JSON} |  awk '{A=gsub(/_/,X,$0)} END {print A}'`
						if [ ${FILEINDEX} == 3 ] ; then
							newfile_json=`basename $JSON | cut -d '_' -f1-2,4`
						elif [ ${FILEINDEX} == 4 ] ; then
							newfile_json=`basename $JSON | cut -d '_' -f1-3,5`
						fi
						mv $JSON $DIR_LOCAL_FMAP/$newfile_json
						oldfile_nifti=`echo $JSON | sed s@'json'@'nii.gz'@g`
						newfile_nifti=`echo $newfile_json | sed s@'json'@'nii.gz'@g`
						mv $oldfile_nifti $DIR_LOCAL_FMAP/$newfile_nifti > /dev/null 2>&1
					fi
				fi
			done
		done
		if [[ ! "$(ls -A $DIR_LOCAL_FMAP)" ]] ; then 
			rm -rf $DIR_LOCAL_FMAP
		fi
	fi

######################################
### Reorganize Directory Structure ###
######################################

	echo "Reorganizing Directory Structure for ${SUBID}"
	chmod -R ug+wrx  ${DIR_PROJECT}/${SUBID}
	if [[ $OPT_RM_DICOMS == TRUE ]] ; then
		rm -rf ${DIR_PROJECT}/${SUBID}/dicoms
	fi
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
