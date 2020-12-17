#!/usr/bin/env bash
#PBS -o /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stdout
#PBS -e /scratch/rjirsara/.submitted_jobs/${PBS_JOBNAME}.stderr
#PBS -l nodes=1:ppn=4,mem=8gb,walltime=24:00:00
#PBS -N dque
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=$1
ATLAS_FILE=$2

module load "freesurfer-5.3.0"

#########
### Extract Freesurfer Data Using Custom Atlas
#########

export SUBJECTS_DIR=${STUDY_DIR}/apps/freesurfer
ATLAS_LABEL=$(basename $ATLAS_FILE | cut -d '.' -f1) ; mkdir -p ${SUBJECTS_DIR}/atlases
if [[ ! -f `echo ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}_*.mgh | cut -d ' ' -f1` ]] ; then
	fslregister \
		--s fsaverage \
		--mov ${ATLAS_FILE} \
		--dof 12 \
		--reg ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}.dat \
	mri_vol2surf \
		--mov ${ATLAS_FILE} \
		--reg ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}.dat \
		--projdist-max 0 1 0.1 \
		--interp nearest \
		--hemi lh \
		--out ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}_lh.mgh
	mri_vol2surf \
		--mov ${ATLAS_FILE} \
		--reg ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}.dat \
		--projdist-max 0 1 0.1 \
		--interp nearest \
		--hemi rh \
		--out ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}_rh.mgh
fi

for SUB_DIR in `ls $SUBJECTS_DIR | grep sub | grep -v base$ | grep -v ERROR$` ; do
	for NEURO in `echo lh.thickness rh.thickness lh.area rh.area lh.volume rh.volume` ; do
		if [[ ! -f `echo ${SUBJECTS_DIR}/${SUB_DIR}/stats/ATL_${ATLAS_LABEL}_${NEURO}.txt` ]] ; then
			HEMI=`echo $NEURO | cut -d '.' -f1`
			mri_surf2surf \
				--s ${SUB_DIR} \
				--trgsubject fsaverage \
				--hemi ${HEMI} \
				--sval ${SUBJECTS_DIR}/${SUB_DIR}/surf/${NEURO} \
				--tval ${SUBJECTS_DIR}/${SUB_DIR}/surf/${NEURO}.fsaverage.mgh
			mri_segstats \
				--seg ${SUBJECTS_DIR}/atlases/${ATLAS_LABEL}_${HEMI}.mgh \
				--in  ${SUBJECTS_DIR}/${SUB_DIR}/surf/${NEURO}.fsaverage.mgh \
				--excludeid 0 \
				--sum ${SUBJECTS_DIR}/${SUB_DIR}/stats/ATL_${ATLAS_LABEL}_${NEURO}.txt
		fi
	done
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
