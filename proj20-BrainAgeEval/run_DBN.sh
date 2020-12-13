#!/usr/bin/env bash
#PBS -l nodes=1:ppn=1,mem=4gb,walltime=4:00:00
#PBS -N MASS_PIPELINE
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=/scratch/rjirsara/study-PDS
PROJECT_DIR=/scratch/rjirsara/proj20-BrainAgeEval
SCRIPTS_DIR=/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval

if [[ ! -d $SCRIPTS_DIR/DeepBrainNet ]] ; then
	git clone https://github.com/RobertJirsaraie/DeepBrainNet.git
fi

#########
### Extract The Brain and Standardize the Intensity Values
#########

for SCAN in `ls $PROJECT_DIR/mass_skullstrip/*_brain.nii.gz` ; do
	SUB=`basename $SCAN | cut -d '_' -f1 | cut -d '-' -f2`
	SES=`basename $SCAN | cut -d '_' -f2 | cut -d '-' -f2`
	if [[ ! -f ${PROJECT_DIR}/brainages_DBN/${SUB}_${SES}/Subject${SUB}${SES}_T1_BrainAligned.nii.gz ]] ; then	
		echo "Extracting the Brain & Standardizing Intensity For SUB: ${SUB} SES: ${SES} "
		mkdir -p $PROJECT_DIR/brainages_DBN/${SUB}_${SES}
		OUTPUT=`echo $PROJECT_DIR/brainages_DBN/${SUB}_${SES}/Subject${SUB}${SES}_T1_BrainAligned.nii.gz`
		flirt -in $SCAN -ref $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152_T1_1mm_brain_LPS_filled.nii.gz -dof 12 -out $OUTPUT
	fi
done
COUNT=`ls ${PROJECT_DIR}/brainages_DBN/*_*/*.nii.gz | wc -l`
fslmerge -t ${PROJECT_DIR}/brainages_DBN/n${COUNT}_skullstripped.nii.gz `ls ${PROJECT_DIR}/brainages_DBN/*_*/*.nii.gz`

#########
### Split the Input Scans Into 80 Sreenshots
#########

for INPUT in `echo ${PROJECT_DIR}/brainages_DBN/*/*.nii.gz` ; do
	if [[ ! -f $(echo `dirname $INPUT`/*79.jpg) ]] ; then
		SUB=`basename $INPUT | cut -d '_' -f1 | cut -d '-' -f2`
		SES=`basename $INPUT | cut -d '_' -f2 | cut -d '-' -f2`
		mkdir -p ${PROJECT_DIR}/${SUB}_${SES}/Test
		python ${SCRIPTS_DIR}/DeepBrainNet/Script/Slicer.py `dirname ${INPUT}`/ ${PROJECT_DIR}/${SUB}_${SES}/
		mv ${PROJECT_DIR}/${SUB}_${SES}/Test/*.jpg `dirname $INPUT` 
		rm -rf ${PROJECT_DIR}/${SUB}_${SES}
	fi
done

#########
### Transfer Input Images to Local Computer Environment Where DeepBrainModel is Configured
#########

#python ${SCRIPTS_DIR}/DeepBrainNet/Script/Model_Test.py \
#	${PROJECT_DIR}/brainages_DBN \
#	${PROJECT_DIR}/brainages_DBN/output.csv \
#	${SCRIPTS_DIR}/DeepBrainNet/CBICA/DBN_model.h5

#LOCAL_DIR=/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval
#rsync -auv rjirsara@login.chpc.wustl.edu:/scratch/rjirsara/proj20-BrainAgeEval/brainages_DBN ${LOCAL_DIR}
#python3 ${LOCAL_DIR}/DeepBrainNet/Script/Model_Test.py ${LOCAL_DIR}/brainages_DBN ${LOCAL_DIR}/brainages_DBN/output.csv ${LOCAL_DIR}/DeepBrainNet/CBICA/DBN_model.h5
#rsync -auv ${LOCAL_DIR}/brainages_DBN/output.csv rjirsara@login.chpc.wustl.edu:/scratch/rjirsara/proj20-BrainAgeEval/brainages_DBN

#########
### Edit Output File To Be Aggregated with the Final Spreadsheet
#########

NSUBS=`cat ${PROJECT_DIR}/brainages_DBN/output.csv | wc -l` ; NSUBS=`echo $(($NSUBS-1))`
cat ${PROJECT_DIR}/brainages_DBN/output.csv | tr '_' ',' | tr '/' ',' | cut -d ',' -f1,2,4 > temp.csv
cat temp.csv | sed s@'ID,Pred'@'sub,ses,brainage_DBN'@g > ${PROJECT_DIR}/brainages_DBN/n${NSUBS}_DeepBrainNet.csv
rm ${PROJECT_DIR}/brainages_DBN/output.csv temp.csv

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
