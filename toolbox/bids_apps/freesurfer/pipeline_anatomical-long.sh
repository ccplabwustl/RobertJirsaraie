#!/usr/bin/env bash
#PBS -o /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stdout
#PBS -e /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stderr
#PBS -l nodes=1:ppn=4,mem=8gb,walltime=24:00:00
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=$1
SUB=$2
SES=$3
module load "freesurfer-5.3.0"

#########
### Process Data Longitudinally
#########

for INPUT in `find $STUDY_DIR/apps/freesurfer -maxdepth 1 | grep "/sub-${SUB}_ses-${SES}" | grep -v ERROR | grep -v long` ; do
	SUB_IDS=`basename $INPUT`	
	export SUBJECTS_DIR=$STUDY_DIR/apps/freesurfer
	BASE_DIR=$STUDY_DIR/apps/freesurfer/sub-${SUB}_base	
	LOG_FILE=$STUDY_DIR/apps/freesurfer/logs/${SUB_IDS}_long.txt
	recon-all -long ${INPUT} ${BASE_DIR} -all > $LOG_FILE
	mv $STUDY_DIR/apps/freesurfer/${SUB_IDS}.long.* $STUDY_DIR/apps/freesurfer/${SUB_IDS}_long
	if [[ -z `cat ${LOG_FILE} | grep 'finished without error' | sed s@' '@'_'@g` ]] ; then
		mv $STUDY_DIR/apps/freesurfer/${SUB_IDS}_long $STUDY_DIR/apps/freesurfer/${SUB_IDS}_long_ERROR
	else
		chmod -R 750 $STUDY_DIR/apps/freesurfer/${SUB_IDS}_long
	fi
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
