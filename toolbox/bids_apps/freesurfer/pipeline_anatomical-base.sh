#!/usr/bin/env bash
#PBS -o /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stdout
#PBS -e /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stderr
#PBS -l nodes=1:ppn=4,mem=8gb,walltime=24:00:00
#PBS -N FSBASE$1
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=$1
SUB=$2
module load "freesurfer-5.3.0"

#########
### Create Base Subject-Specific Templates
#########

mkdir -p $STUDY_DIR/apps/freesurfer/logs 
export SUBJECTS_DIR=$STUDY_DIR/apps/freesurfer
LOG_FILE=$STUDY_DIR/apps/freesurfer/logs/sub-${SUB}_base.txt
CROSS=`echo $STUDY_DIR/apps/freesurfer/sub-${SUB}_ses-* | tr ' ' '\n' | grep -v ERROR`
ROOTPATH=`echo $CROSS | cut -d ' ' -f1 | cut -d '/' -f2`
INPUTS=`echo $CROSS | sed s@/"${ROOTPATH}"@"-tp /${ROOTPATH}"@g`
recon-all -base $STUDY_DIR/apps/freesurfer/sub-${SUB}_base ${INPUTS} -all -qcache > $LOG_FILE
if [[ -z `cat ${LOG_FILE} | grep 'finished without error' | sed s@' '@'_'@g` ]] ; then
	mv $STUDY_DIR/apps/freesurfer/sub-${SUB}_base $STUDY_DIR/apps/freesurfer/sub-${SUB}_base_ERROR
else
	chmod -R 750 $STUDY_DIR/apps/freesurfer/sub-${SUB}_base
fi

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
