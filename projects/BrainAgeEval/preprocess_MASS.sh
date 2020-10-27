#!/bin/bash
#PBS -o /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stdout
#PBS -e /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stderr
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=12:00:00
#PBS -N MASS_PIPELINE
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=/scratch/rjirsara/study-PDS
PROJECT_DIR=/scratch/rjirsara/projects/BrainAgeEval/DeepBrainNet
SCRIPTS_DIR=/scratch/rjirsara/RobertJirsaraie/projects/BrainAgeEval
module purge ; module load "fsl-5.0.8" "dramms-1.5.1" "AFNI-17.0.9"

##########

mkdir -p $PROJECT_DIR/mass_N4correct
cat $SCRIPTS_DIR/mass_1.1.1/lib/mass-jobsubmit-PDS > $PROJECT_DIR/mass_N4correct/command.sh
echo 'for RAW in `find $STUDY_DIR/bids/sub-*/ses-*/anat/*_T1w.nii.gz | grep -v run` ; do' >> $PROJECT_DIR/mass_N4correct/command.sh
echo '	if [[ ! -f $(echo $PROJECT_DIR/mass_N4correct/`basename $RAW`) ]] ; then' >> $PROJECT_DIR/mass_N4correct/command.sh
echo '		module purge ; module load "ANTs"' >> $PROJECT_DIR/mass_N4correct/command.sh
echo '		N4BiasFieldCorrection -i $RAW -o $PROJECT_DIR/mass_N4correct/`basename $RAW`' >> $PROJECT_DIR/mass_N4correct/command.sh
echo '	fi ' >> $PROJECT_DIR/mass_N4correct/command.sh
echo 'done ' >> $PROJECT_DIR/mass_N4correct/command.sh
echo '' >> $PROJECT_DIR/mass_N4correct/command.sh
chmod -R ug+wrx $PROJECT_DIR/mass_N4correct
qsub -N "N4CORRECT" $PROJECT_DIR/mass_N4correct/command.sh

####################################################################
### Create Study Specific Templates For Each of The Two Scanners ###
####################################################################

COUNT=35
ROOT_DIR=$PROJECT_DIR/mass_kclusters${COUNT} ; mkdir -p $ROOT_DIR
cat $SCRIPTS_DIR/mass_1.1.1/lib/mass-jobsubmit-PDS > $ROOT_DIR/command.sh
find $PROJECT_DIR/mass_N4correct/*.nii.gz > $ROOT_DIR/n837_kcohort.lst
EXP=`cat $ROOT_DIR/n837_kcohort.lst | head -n1` 
3dresample -dxyz 1.0 1.0 1.0 -input ${EXP} -prefix test.nii.gz ; mv test.nii.gz ${EXP}
echo "EXP=$EXP" >> $ROOT_DIR/command.sh
echo "COUNT=$COUNT" >> $ROOT_DIR/command.sh
echo 'for FILE in `cat $ROOT_DIR/n837_kcohort.lst` ; do' >> $ROOT_DIR/command.sh
echo 	'flirt -in $FILE -ref ${EXP} -out $FILE -applyxfm' >> $ROOT_DIR/command.sh
echo 'done' >> $ROOT_DIR/command.sh
echo "" >> $ROOT_DIR/command.sh
echo '$SCRIPTS_DIR/bin/mass-chooseTemplates \
	-list $ROOT_DIR/n837_kcohort.lst \
	-dest $ROOT_DIR/ \
	-tmp $ROOT_DIR/ \
	-clust ${COUNT} \
	-MT 8 -v' >> $ROOT_DIR/command.sh
echo "" >> $ROOT_DIR/command.sh
echo "" >> $ROOT_DIR/command.sh
echo "chmod ug+wrx $ROOT_DIR/command.sh" >> $ROOT_DIR/command.sh
echo "rm $ROOT_DIR/ChosenTemplates_*_metrics.txt" >> $ROOT_DIR/command.sh
echo "mv $ROOT_DIR/ChosenTemplates_*.txt $ROOT_DIR/n${COUNT}_kclusters.lst" >> $ROOT_DIR/command.sh
qsub -N "MASSCLUST" $ROOT_DIR/command.sh


##########

for TEMPLATE in `cat $PROJECT_DIR/mass_kclusters35/n35_kclusters.lst` ; do
	$SCRIPTS_DIR/mass_1.1.1/bin/mass \
		-in $PROJECT_DIR/mass_N4correct/`basename $TEMPLATE` \
		-dest $PROJECT_DIR/mass_templates/ \
		-regs 15 \
		-agg 50 \
		-int 0 \
		-MT 6 
done

##########

INDEX=0
for SCAN in `ls $PROJECT_DIR/mass_templates/sub-*_ses-*_T1w.nii.gz` ; do
	INDEX=$(($INDEX+1))
	mv $SCAN $(dirname $SCAN)/Template${INDEX}.nii.gz
done
INDEX=0
for SCAN in `ls $PROJECT_DIR/mass_templates/sub-*_ses-*_T1w_brainmask.nii.gz` ; do
	INDEX=$(($INDEX+1))
	mv $SCAN $(dirname $SCAN)/Template${INDEX}_str_cbq.nii.gz
done
INDEX=0
for SCAN in `ls $PROJECT_DIR/mass_templates/sub-*_ses-*_T1w_brain.nii.gz` ; do
	INDEX=$(($INDEX+1))
	mv $SCAN $(dirname $SCAN)/Template${INDEX}_brain.nii.gz
done
rm `ls $PROJECT_DIR/mass_templates/* | grep -v /Template`

######

mkdir -p $PROJECT_DIR/mass_skullstrip
for RAW in `cat $ROOT_DIR/n837_kcohort.lst | head -n3` ; do
	$SCRIPTS_DIR/mass_1.1.1/bin/mass \
		-in $RAW \
		-ref $PROJECT_DIR/mass_templates \
		-dest $PROJECT_DIR/mass_skullstrip \
		-regs 20 \
		-agg 50 \
		-int 0 \
		-MT 6 
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######

