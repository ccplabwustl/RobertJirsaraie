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
PROJECT_DIR=/scratch/rjirsara/proj20-BrainAgeEval
SCRIPTS_DIR=/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval
module purge ; module load "fsl-5.0.8" "dramms-1.5.1" "AFNI-17.0.9"

#########
### Apply Bias Correction to Raw Scans
#########

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
qsub -N "N4CORRECT" $PROJECT_DIR/mass_N4correct/command.sh ; mv $PROJECT_DIR/mass_N4correct/command.sh $PROJECT_DIR/mass_N4correct/commandcor.sh 

#########
### Standardized Image Resolution, Number of Voxels, and Normalize to MNI Space
#########

TEMPLATE=$SCRIPTS_DIR/mass_1.1.1/data/
cp ~/.cache/templateflow/tpl-MNI152Lin/tpl-MNI152Lin_res-01_T1w.nii.gz $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_T1w.nii.gz
3dresample -input $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_T1w.nii.gz \
	-prefix $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_TEMP1_T1w.nii.gz \
	-dxyz 1 1 1
flirt -in $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_TEMP1_T1w.nii.gz \
	-ref `ls $PROJECT_DIR/mass_N4correct/*.nii.gz | head -n1` \
	-out $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_TEMP2_T1w.nii.gz \
	-applyxfm
mv $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_TEMP2_T1w.nii.gz $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_T1w.nii.gz
rm $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_TEMP1_T1w.nii.gz
cat $SCRIPTS_DIR/mass_1.1.1/lib/mass-jobsubmit-PDS > $PROJECT_DIR/mass_N4correct/commandreg.sh
echo 'for SCAN in `ls $PROJECT_DIR/mass_N4correct/*.nii.gz` ; do' >> $PROJECT_DIR/mass_N4correct/commandreg.sh
echo '	flirt -in ${SCAN} -ref $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_T1w.nii.gz -dof 12 -out test.nii.gz' >> $PROJECT_DIR/mass_N4correct/commandreg.sh
echo '	mv test.nii.gz ${SCAN} ' >> $PROJECT_DIR/mass_N4correct/commandreg.sh
echo 'done ' >> $PROJECT_DIR/mass_N4correct/commandreg.sh
echo '' >> $PROJECT_DIR/mass_N4correct/commandreg.sh
chmod -R ug+wrx $PROJECT_DIR/mass_N4correct
qsub -N "FLIRTREG" $PROJECT_DIR/mass_N4correct/commandreg.sh

#########
### Run K-means to Select Study-Specific Templates
#########

COUNT=6
mkdir -p $PROJECT_DIR/mass_kclusters
cat $SCRIPTS_DIR/mass_1.1.1/lib/mass-jobsubmit-PDS > $PROJECT_DIR/mass_kclusters/command.sh
find $PROJECT_DIR/mass_N4correct/*.nii.gz > $PROJECT_DIR/mass_kclusters/n834_kcohort.lst
echo "COUNT=$COUNT" >> $PROJECT_DIR/mass_kclusters/command.sh
echo "$SCRIPTS_DIR/mass_1.1.1/bin/mass-chooseTemplates \
	-list $PROJECT_DIR/mass_kclusters/n834_kcohort.lst \
	-dest $PROJECT_DIR/mass_kclusters/ \
	-tmp $PROJECT_DIR/mass_kclusters/ \
	-clust ${COUNT} \
	-MT 8 -v" >> $PROJECT_DIR/mass_kclusters/command.sh
echo "" >> $PROJECT_DIR/mass_kclusters/command.sh
echo "" >> $PROJECT_DIR/mass_kclusters/command.sh
echo "chmod ug+wrx $PROJECT_DIR/mass_kclusters/command.sh" >> $PROJECT_DIR/mass_kclusters/command.sh
echo "rm $PROJECT_DIR/mass_kclusters/ChosenTemplates_*_metrics.txt" >> $PROJECT_DIR/mass_kclusters/command.sh
echo "mv $PROJECT_DIR/mass_kclusters/ChosenTemplates_*.txt $PROJECT_DIR/mass_kclusters/n${COUNT}_kclusters.lst" >> $PROJECT_DIR/mass_kclusters/command.sh
qsub -N "MASSCLUST" $PROJECT_DIR/mass_kclusters/command.sh

#########
### Register the PNC Templates to Create Study-Specific Templates 
#########

COUNT=6
mkdir -p $PROJECT_DIR/mass_templates
for TEMPLATE in `cat $PROJECT_DIR/mass_kclusters/n${COUNT}_kclusters.lst` ; do
	$SCRIPTS_DIR/mass_1.1.1/bin/mass \
	-in $PROJECT_DIR/mass_N4correct/`basename $TEMPLATE` \
	-dest $PROJECT_DIR/mass_templates/ \
	-regs 6 \
	-agg 50 \
	-int 2 \
	-MT 6 
done

#########
### Relabel the Study Specific Templates To Prepare for SkullStripping
#########

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

#########
### Execute SkullStripping for all Scans based on Study Specific Templates
#########

mkdir -p $PROJECT_DIR/mass_skullstrip
for RAW in `cat $PROJECT_DIR/mass_kclusters/n834_kcohort.lst` ; do
	SUBIDS=`basename $RAW | cut -d '_' -f1,2`
	if [[ ! -f $PROJECT_DIR/mass_skullstrip/${SUBIDS}_T1w_brainmask.nii.gz ]] ; then
	$SCRIPTS_DIR/mass_1.1.1/bin/mass \
		-in $RAW \
		-ref $PROJECT_DIR/mass_templates \
		-dest $PROJECT_DIR/mass_skullstrip \
		-regs 6 \
		-agg 50 \
		-int 0 \
		-MT 6 
	fi
done

#########
### Remove Incomplete Processing and Concatenate All Brain Images for QAing
#########

for INCOMPLETE in `ls -l ${PROJECT_DIR}/mass_skullstrip/*affine* | awk {'print $9'} | cut -d '_' -f1,2 | uniq` ; do
	rm -rf ${INCOMPLETE}*
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
