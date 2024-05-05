#!/usr/bin/env bash
#PBS -o /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stdout
#PBS -e /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stderr
#PBS -l nodes=1:ppn=8,mem=8gb,walltime=24:00:00
#PBS -N replicate_devhcp
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=/scratch/rjirsara/study-PDS
PROJECT_DIR=/scratch/rjirsara/proj20-BrainAgeEval
HCP_DIR=/NRG-data/NRG/intradb/archive/CCF_HCD_STG/arc001
SCRIPTS_DIR=/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval
module purge ; module load "fsl-5.0.8" "dramms-1.5.1" "AFNI-17.0.9" "ANTs" "freesurfer" "gcc-7.5.0"
export PATH=$PATH:${SCRIPTS_DIR}/hcp_workbench/bin_linux64 ; export SUBJECTS_DIR=$STUDY_DIR/apps/freesurfer

#####
### Extract the Euler Number For Each Subject
#####

mkdir -p $PROJECT_DIR/replicate_devhcp/EulerNumber
for SUB in `cat ${PROJECT_DIR}/replicate_devhcp/n789_RAW.csv  | awk -F "\"*,\"*" '{print $1}' | grep -v sub` ; do
	if [[ ! -f $PROJECT_DIR/replicate_devhcp/EulerNumber/${SUB}.txt ]] ; then
		SUBDIR=$HCP_DIR/$SUB/RESOURCES/Structural_preproc/$SUB/T1w/$SUB
		RIGHT=`script -c "mris_euler_number ${SUBDIR}/surf/rh.orig.nofix" | grep ">" | cut -d ">" -f 1 | cut -d "=" -f 4 | cut -d " " -f 2`
		LEFT=`script -c "mris_euler_number ${SUBDIR}/surf/lh.orig.nofix" | grep ">" | cut -d ">" -f 1 | cut -d "=" -f 4 | cut -d " " -f 2`
		SUM=`echo $(( ${LEFT} + ${RIGHT} ))` ; AVG=`echo $(($SUM/2))` ; echo $AVG > $PROJECT_DIR/replicate_devhcp/EulerNumber/${SUB}.txt
		rm typescript ; echo "Finishing Extraction for: `basename $SUBDIR`"
	fi
done

#####
### Extract Regional Features from Glasser 2016 Atlas
#####

mkdir -p $PROJECT_DIR/replicate_devhcp/GlasserLabel
ANNOT_LH=`ls $PROJECT_DIR/freesurf_extract/lh.HCP-MMP1.annot` ; ANNOT_RH=`ls $PROJECT_DIR/freesurf_extract/rh.HCP-MMP1.annot` 
for SUB in `cat ${PROJECT_DIR}/replicate_devhcp/n789_RAW.csv  | awk -F "\"*,\"*" '{print $1}' | grep -v sub` ; do
	if [[ ! -f `echo $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_lh.annot` ]] ; then
		cp -r $HCP_DIR/$SUB/RESOURCES/Structural_preproc/$SUB/T1w/$SUB ${SUBJECTS_DIR}
		mri_surf2surf \
			--srcsubject fsaverage \
			--trgsubject $SUB \
			--hemi lh \
			--sval-annot ${ANNOT_LH} \
			--tval $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_lh.annot
		mri_surf2surf \
			--srcsubject fsaverage \
			--trgsubject $SUB \
			--hemi rh \
			--sval-annot ${ANNOT_RH} \
			--tval $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_rh.annot
		mris_anatomical_stats \
			-mgz -cortex $SUBJECTS_DIR/${SUB}/label/lh.cortex.label \
			-f $SUBJECTS_DIR/${SUB}/stats/lh.HCP-MMP1.stats -b \
			-a $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_lh.annot ${SUB} lh white
		cp $SUBJECTS_DIR/${SUB}/stats/lh.HCP-MMP1.stats $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_corticallh.stats
		mris_anatomical_stats \
			-mgz -cortex $SUBJECTS_DIR/${SUB}/label/rh.cortex.label \
			-f $SUBJECTS_DIR/${SUB}/stats/rh.HCP-MMP1.stats -b \
			-a $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_rh.annot ${SUB} rh white
		cp $SUBJECTS_DIR/${SUB}/stats/rh.HCP-MMP1.stats $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_corticalrh.stats
		cp $SUBJECTS_DIR/${SUB}/stats/aseg.stats $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_subcort.stats
		rm -rf $SUBJECTS_DIR/${SUB}
	fi
done

#####
### Aggregate Regional Features Into Master Spreadsheet
#####

for SUB in `ls $PROJECT_DIR/replicate_devhcp/GlasserLabel/*lh.stats | cut -d '/' -f7 | cut -d '_' -f1,2,3` ; do
	if [[ ! -f $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}.csv ]] ; then
		for HEMI in `echo lh rh` ; do
			cat $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_cortical${HEMI}.stats | head -n241 | tail -n 180 > ${SUB}_${HEMI}.csv
			cat ${SUB}_${HEMI}.csv | awk {'print $1,$3'} | tr ' ' ',' | sed s@'ROI'@'ROI_area'@g >> ${SUB}.csv
			cat ${SUB}_${HEMI}.csv | awk {'print $1,$4'} | tr ' ' ',' | sed s@'ROI'@'ROI_volume'@g >> ${SUB}.csv
			cat ${SUB}_${HEMI}.csv | awk {'print $1,$5'} | tr ' ' ',' | sed s@'ROI'@'ROI_thickness'@g >> ${SUB}.csv
		done
		FILE=`ls $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}_subcort.stats`
		cat ${SUB}.csv | sed s@"^R_"@"rh_R_"@g | sed s@"^L_"@"lh_L_"@g > ${SUB}_TEMP.csv ; mv ${SUB}_TEMP.csv ${SUB}.csv
		cat $FILE | grep -v ^#  | awk {'print $5,$4'} | tr ' ' ',' | sed s@'-'@'.'@g >> ${SUB}.csv
		cat ${SUB}.csv | grep -v hypointensities  | grep -v vessel | grep -v Optic | grep -v VentralDC | grep -v choroid > ${SUB}_TEMP.csv
		cat ${SUB}_TEMP.csv | sed s@'3rd'@'X3rd'@g | sed s@'4th'@'X4th'@g | grep -v '5th.Ventricle' | grep -v 'CSF' > ${SUB}_TEMP1.csv
		echo lhCortexVol,`cat $FILE | grep lhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo rhCortexVol,`cat $FILE | grep rhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo lhCorticalWhiteMatterVol,`cat $FILE | grep lhCerebralWhiteMatterVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo rhCorticalWhiteMatterVol,`cat $FILE | grep rhCerebralWhiteMatterVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo SubCortGrayVol,`cat $FILE | grep SubCortGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo TotalGrayVol,`cat $FILE | grep TotalGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo SupraTentorialVol,`cat $FILE | grep SupraTentorialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		echo EstimatedTotalIntraCranialVol,`cat $FILE | grep EstimatedTotalIntraCranialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUB}_TEMP1.csv
		cat ${SUB}_TEMP1.csv | sed s@'-'@'.'@g > ${SUB}_TEMP2.csv ; rm `ls ${SUB}*.csv | grep -v _TEMP2.csv`
		mv ${SUB}_TEMP2.csv $PROJECT_DIR/replicate_devhcp/GlasserLabel/${SUB}.csv
	fi
done

for SUBFILE in `ls $PROJECT_DIR/replicate_devhcp/GlasserLabel/*.csv` ; do
	if [[ `cat $SUBFILE | wc -l` != 1118 ]] ; then
		echo "Incomplete Dataset Removed: `basename $SUBFILE` - `cat $SUBFILE | wc -l`" ; rm $SUBFILE
	fi
done

#########
### Apply Bias Correction and Normalize all Raw Scans
#########

mkdir -p $PROJECT_DIR/replicate_devhcp/CopiedT1w $PROJECT_DIR/replicate_devhcp/MASS
for RAW in `cat ${PROJECT_DIR}/replicate_devhcp/n789_RAW.csv | tr ',' '\n'  | grep HCD` ; do
	if [[ -f /NRG-data/NRG/intradb/archive/CCF_HCD_STG/arc001/${RAW}/RESOURCES/Structural_preproc/${RAW}/T1w/T1w.nii.gz ]] ; then
		cp /NRG-data/NRG/intradb/archive/CCF_HCD_STG/arc001/${RAW}/RESOURCES/Structural_preproc/${RAW}/T1w/T1w.nii.gz ${PROJECT_DIR}/replicate_devhcp/CopiedT1w/${RAW}_T1w.nii.gz
		N4BiasFieldCorrection \
			-i ${PROJECT_DIR}/replicate_devhcp/CopiedT1w/${RAW}_T1w.nii.gz \
			-o ${PROJECT_DIR}/replicate_devhcp/CopiedT1w/${RAW}_BiasCor.nii.gz
		flirt \
			-in ${PROJECT_DIR}/replicate_devhcp/CopiedT1w/${RAW}_BiasCor.nii.gz \
			-ref $SCRIPTS_DIR/mass_1.1.1/data/tpl-MNI152Lin_res-01_T1w.nii.gz \
			-dof 12 \
			-out ${PROJECT_DIR}/replicate_devhcp/CopiedT1w/${RAW}_NormCor.nii.gz
	fi
done

#########
### Run K-means to Select Study-Specific Templates
#########

COUNT=6 ; mkdir -p $PROJECT_DIR/replicate_devhcp/KSelect
cat $SCRIPTS_DIR/mass_1.1.1/lib/mass-jobsubmit-PDS > $PROJECT_DIR/replicate_devhcp/KClusters/TemplateSelect.sh
find $PROJECT_DIR/replicate_devhcp/CopiedT1w/*_NormCor.nii.gz > $PROJECT_DIR/replicate_devhcp/KClusters/n789_kcohort.lst
echo "$SCRIPTS_DIR/mass_1.1.1/bin/mass-chooseTemplates \
	-list $PROJECT_DIR/replicate_devhcp/KClusters/n789_kcohort.lst \
	-dest $PROJECT_DIR/replicate_devhcp/KClusters/ \
	-tmp $PROJECT_DIR/replicate_devhcp/KClusters \
	-clust ${COUNT} \
	-MT 8 -v" >> $PROJECT_DIR/replicate_devhcp/KClusters/TemplateSelect.sh
qsub -N "MASSCLUST" $PROJECT_DIR/replicate_devhcp/KClusters/TemplateSelect.sh
rm $PROJECT_DIR/replicate_devhcp/KClusters/ChosenTemplates_*_metrics.txt
mv $PROJECT_DIR/replicate_devhcp/KClusters/ChosenTemplates_*.txt $PROJECT_DIR/replicate_devhcp/KClusters/n${COUNT}_kclusters.lst

#########
### Register the PNC Templates to Create Study-Specific Templates 
#########

COUNT=6 ; mkdir -p $PROJECT_DIR/replicate_devhcp/KTemplates
for TEMPLATE in `cat $PROJECT_DIR/replicate_devhcp/KSelect/n${COUNT}_kclusters.lst` ; do
	$SCRIPTS_DIR/mass_1.1.1/bin/mass \
		-in $TEMPLATE \
		-dest $PROJECT_DIR/replicate_devhcp/KTemplates \
		-regs 6 \
		-agg 50 \
		-int 2 \
		-MT 6 
done

#########
### Execute SkullStripping for all Scans based on Study Specific Templates
#########

mkdir -p $PROJECT_DIR/replicate_devhcp/MASS
for RAW in `cat $PROJECT_DIR/replicate_devhcp/KSelect/n789_kcohort.lst` ; do
	SUBIDS=`basename $RAW  | sed s@.nii.gz@@g`
	if [[ ! -f $PROJECT_DIR/replicate_devhcp/MASS/${SUBIDS}_brain.nii.gz ]] ; then
		$SCRIPTS_DIR/mass_1.1.1/bin/mass \
			-in $RAW \
			-ref $PROJECT_DIR/mass_templates \
			-dest $PROJECT_DIR/replicate_devhcp/MASS \
			-regs 6 \
			-agg 50 \
			-int 0 \
			-MT 6 
	fi
done

#########
### Split the Input Scans Into 80 Sreenshots
#########

mkdir -p $PROJECT_DIR/replicate_devhcp/Test
python $SCRIPTS_DIR/DeepBrainNet/Script/Slicer.py \
	$PROJECT_DIR/replicate_devhcp/MASS/ \
	$PROJECT_DIR/replicate_devhcp
 mv $PROJECT_DIR/replicate_devhcp/Test $PROJECT_DIR/replicate_devhcp/DBNInput

#########
### Execute SkullStripping for all Scans based on Study Specific Templates
#########

python ${SCRIPTS_DIR}/DeepBrainNet/Script/Model_Test.py \
	${PROJECT_DIR}/replicate_devhcp/DBNInput \
	${PROJECT_DIR}/replicate_devhcp/DBNInput/DBN.csv \
	${SCRIPTS_DIR}/DeepBrainNet/Models/DBN_model.h5

python ${SCRIPTS_DIR}/DeepBrainNet/Script/Model_Test.py \
	${PROJECT_DIR}/replicate_devhcp/DBNInput \
	${PROJECT_DIR}/replicate_devhcp/DBNInput/rDBN.csv \
	${SCRIPTS_DIR}/DeepBrainNet/Models/DBN_Pediatric.h5

python3 /Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/DeepBrainNet/Script/Model_Test.py \
/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/replicate_devhcp/DBNInput \
/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/replicate_devhcp/DBN.csv \
/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/DeepBrainNet/Models/DBN_Pediatric.h5

#########
### Edit Output File To Be Aggregated with the Final Spreadsheet
#########

NSUBS=`cat ${PROJECT_DIR}/analysis/brainages_DBN/output.csv | wc -l` ; NSUBS=`echo $(($NSUBS-1))`
cat ${PROJECT_DIR}/analysis/brainages_DBN/output.csv | tr '_' ',' | tr '/' ',' | cut -d ',' -f1,2,4 > temp.csv
cat temp.csv | sed s@'ID,Pred'@'sub,ses,brainage_DBN'@g > ${PROJECT_DIR}/analysis/brainages_DBN/n${NSUBS}_DeepBrainNet.csv
rm ${PROJECT_DIR}/analysis/brainages_DBN/output.csv temp.csv

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####          ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡          ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
