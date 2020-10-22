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
module load "fsl-5.0.8" "dramms-1.5.1" "AFNI-17.0.9"

####################################################################
### Create Study Specific Templates For Each of The Two Scanners ###
####################################################################

COUNT=15
cat $SCRIPTS_DIR/mass_1.1.1/lib/mass-jobsubmit-PDS > $PROJECT_DIR/command_kmeans.sh
find $STUDY_DIR/bids/sub-*/ses-[1..2..3]/anat/*_T1w.nii.gz | \
	grep -v run |  \
	/bin/sort -R | \
	head -n554 > $PROJECT_DIR/n554_kcohort123.lst
echo "$SCRIPTS_DIR/mass_1.1.1/bin/mass-chooseTemplates \
	-list $PROJECT_DIR/n554_kcohort123.lst \
	-dest $PROJECT_DIR/ \
	-tmp $PROJECT_DIR/ \
	-clust ${COUNT} \
	-MT 8 -v" >> $PROJECT_DIR/command_kmeans.sh
chmod ug+wrx $PROJECT_DIR/command_kmeans.sh
qsub -N "MASSCLUST" $PROJECT_DIR/command_kmeans.sh
rm $PROJECT_DIR/ChosenTemplates_*_metrics.txt ; mv $PROJECT_DIR/ChosenTemplates_*.txt $PROJECT_DIR/n${COUNT}_kclusters123.lst

##########

for TEMPLATE in `cat $PROJECT_DIR/n${COUNT}_kclusters123.lst` ; do
	TEMPLATE=/scratch/rjirsara/study-PDS/bids/sub-606/ses-2/anat/sub-606_ses-2_T1w.nii.gz
	$SCRIPTS_DIR/mass_1.1.1/bin/mass -in ${TEMPLATE} -dest $PROJECT_DIR/mass_template123/run3/ -NOQ -MT 6 -regs 15 -int 2

-regs	  < int   >	No. of templates
 -int 2
done


TEMPLATE=/scratch/rjirsara/study-PDS/bids/sub-404/ses-2/anat/sub-404_ses-2_T1w.nii.gz

$SCRIPTS_DIR/bin/mass -in ${TEMPLATE} -dest $PROJECT_DIR/mass_template123 -NOQ -MT 6







find $STUDY_DIR//bids/sub-*/ses-[4..5]/anat/*_T1w.nii.gz | \
	grep -v run |  \
	/bin/sort -R | \
	head -n25 > $PROJECT_DIR/mass_template45/choosetemplate.lst
COMMAND_FILE=`echo $PROJECT_DIR/mass_template123/command.sh`
cat $SCRIPTS_DIR//lib/mass-jobsubmit > $COMMAND_FILE
echo "$SCRIPTS_DIR/bin/mass-chooseTemplates -list $PROJECT_DIR/mass_template123/choosetemplate.lst -clust 15 -tmp $PROJECT_DIR/mass_template123 -MT 8 -v" >> $COMMAND_FILE
chmod ug+wrx $COMMAND_FILE ; qsub $COMMAND_FILE


########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######

