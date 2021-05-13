#!/usr/bin/env bash
#PBS -o /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stdout
#PBS -e /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stderr
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=48:00:00
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=$1
SUB=$2
SES=$3
module purge ; module load "freesurfer-5.3.0"

#########
### Process Cross-sectional Data To Get Base Subject-Specific Templates
#########

for INPUT in `find $STUDY_DIR/bids | grep T1w.nii.gz | grep "/sub-${SUB}_ses-${SES}"` ; do
	SUB_IDS=`basename $INPUT | sed s@'_T1w.nii.gz'@''@g`
	export SUBJECTS_DIR=$STUDY_DIR/apps/freesurfer
	MGZ_FILE=$STUDY_DIR/apps/freesurfer/logs/${SUB_IDS}_T1w.mgz
	LOG_FILE=$STUDY_DIR/apps/freesurfer/logs/${SUB_IDS}_cross.txt
	if [ ! -f ${MGZ_FILE} ] ; then
		mkdir -p $STUDY_DIR/apps/freesurfer/logs
		mri_convert $INPUT ${MGZ_FILE} > ${LOG_FILE} 
		rm -rf $STUDY_DIR/apps/freesurfer/${SUB_IDS} ; chmod 750 ${MGZ_FILE}
	fi
	recon-all -i ${MGZ_FILE} -s ${SUB_IDS} -all -qcache >> ${LOG_FILE}
	if [[ -z `cat ${LOG_FILE} | grep 'finished without error' | sed s@' '@'_'@g` ]] ; then
		mv $STUDY_DIR/apps/freesurfer/$SUB_IDS $STUDY_DIR/apps/freesurfer/${SUB_IDS}_ERROR
	else
		chmod -R 750 $STUDY_DIR/apps/freesurfer/${SUB_IDS}
	fi
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
