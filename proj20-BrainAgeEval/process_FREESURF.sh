#!/usr/bin/env bash
#PBS -l nodes=1:ppn=8,mem=8gb,walltime=24:00:00
#PBS -N FREESURF_PIPE
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=/scratch/rjirsara/study-PDS
PROJECT_DIR=/scratch/rjirsara/proj20-BrainAgeEval
SCRIPTS_DIR=/scratch/rjirsara/RobertJirsaraie/proj20-BrainAgeEval
export PATH=$PATH:${SCRIPTS_DIR}/hcp_workbench/bin_linux64
export SUBJECTS_DIR=$STUDY_DIR/apps/freesurfer

#####
### Download the correct data from BALSA:
### https://balsa.wustl.edu/WN56
#####

#scp ~/Desktop/Glasser_et_al_2016_HCP_MMP1.0_qN_RVVG.zip rjirsara@login.chpc.wustl.edu:${PROJECT_DIR}/freesurf_extract/
unzip ${PROJECT_DIR}/freesurf_extract/Glasser_et_al_2016_HCP_MMP1.0_qN_RVVG.zip 
rm ${PROJECT_DIR}/freesurf_extract/Glasser_et_al_2016_HCP_MMP1.0_qN_RVVG.zip
mv ${PROJECT_DIR}/freesurf_extract/Glasser_et_al_2016_HCP_MMP1.0_qN_RVVG ${PROJECT_DIR}/freesurf_extract/Glasser

#####
### Download the correct spheres from the HCP github page:
### https://www.mail-archive.com/hcp-users%40humanconnectome.org/msg02890.html
#####

git clone git@github.com:Washington-University/HCPpipelines.git
mv HCPpipelines ${PROJECT_DIR}/freesurf_extract/

#####
### Download and Load Connectome WorkBench Software Package
#####

wget https://www.humanconnectome.org/storage/app/media/workbench/workbench-linux64-v1.4.2.zip
unzip workbench-linux64-v1.4.2.zip ; rm workbench-linux64-v1.4.2.zip ; mv workbench ${SCRIPTS_DIR}/hcp_workbench

#####
### Use wb_command on Local Machine Because of Incompatibility with CHPC
#####

<<RUN_LOCALLY
LOCAL_DIR=/Users/Jirsaraie/Desktop/Research/proj20-BrainAgeEval/freesurf_extract
scp -r rjirsara@login.chpc.wustl.edu:${STUDY_DIR}/apps/freesurfer/fsaverage ${LOCAL_DIR}
export SUBJECTS_DIR=${LOCAL_DIR} ; export PATH=$PATH:${LOCAL_DIR}/workbench/bin_macosx64

NII=${LOCAL_DIR}/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii
GII=${LOCAL_DIR}/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.label.gii
HCP1=${LOCAL_DIR}/HCPpipelines/global/templates/standard_mesh_atlases/L.sphere.32k_fs_LR.surf.gii
HCP2=${LOCAL_DIR}/HCPpipelines/global/templates/standard_mesh_atlases/fs_L/fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii
wb_command -cifti-separate ${NII} COLUMN -label CORTEX_LEFT ${GII}
wb_command -label-resample ${GII} ${HCP1} ${HCP2} BARYCENTRIC left.fsaverage164.label.gii

NII=${LOCAL_DIR}/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii
GII=${LOCAL_DIR}/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.label.gii
HPC1=${LOCAL_DIR}/HCPpipelines/global/templates/standard_mesh_atlases/R.sphere.32k_fs_LR.surf.gii
HPC2=${LOCAL_DIR}/HCPpipelines/global/templates/standard_mesh_atlases/fs_R/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii
wb_command -cifti-separate ${NII} COLUMN -label CORTEX_RIGHT ${GII}
wb_command -label-resample ${GII} ${HCP1} ${HCP2} BARYCENTRIC right.fsaverage164.label.gii

scp -r ${LOCAL_DIR}/*.fsaverage164.label.gii rjirsara@login.chpc.wustl.edu:${PROJECT_DIR}/freesurf_extract
RUN_LOCALLY

#####
### Convert GII File into Annotation Format
#####

GII=$PROJECT_DIR/freesurf_extract/left.fsaverage164.label.gii
HCP=`find $PROJECT_DIR/freesurf_extract/HCPpipelines -iname fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii`
mris_convert --annot ${GII} ${HCP} $PROJECT_DIR/freesurf_extract/lh.HCP-MMP1.annot

GII=$PROJECT_DIR/freesurf_extract/right.fsaverage164.label.gii
HCP=`find $PROJECT_DIR/freesurf_extract/HCPpipelines -iname fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii`
mris_convert --annot ${GII} ${HCP} $PROJECT_DIR/freesurf_extract/rh.HCP-MMP1.annot

#####
### Extract Regional Features from Glasser 2016 Atlas
#####

ANNOT_LH=`ls $PROJECT_DIR/freesurf_extract/lh.HCP-MMP1.annot`
ANNOT_RH=`ls $PROJECT_DIR/freesurf_extract/rh.HCP-MMP1.annot` 
for SUBNAME in `ls $SUBJECTS_DIR | grep sub | grep -v ERROR | grep -v base` ; do
	if [[ ! -f `echo $SUBJECTS_DIR/${SUBNAME}/stats/lh.HCP-MMP1.stats` ]] ; then
		mri_surf2surf \	
			--srcsubject fsaverage \
			--trgsubject ${SUBNAME} \
			--hemi lh \
			--sval-annot ${ANNOT_LH} \
			--tval $SUBJECTS_DIR/${SUBNAME}/label/lh.HCP-MMP1.annot
		mri_surf2surf \
			--srcsubject fsaverage \
			--trgsubject ${SUBNAME} \
			--hemi rh \
			--sval-annot ${ANNOT_RH} \
			--tval $SUBJECTS_DIR/${SUBNAME}/label/rh.HCP-MMP1.annot
		mris_anatomical_stats \
			-mgz -cortex $SUBJECTS_DIR/${SUBNAME}/label/lh.cortex.label \
			-f $SUBJECTS_DIR/${SUBNAME}/stats/lh.HCP-MMP1.stats -b \
			-a $SUBJECTS_DIR/${SUBNAME}/label/lh.HCP-MMP1.annot ${SUBNAME} lh white
		mris_anatomical_stats \
			-mgz -cortex $SUBJECTS_DIR/${SUBNAME}/label/rh.cortex.label \
			-f $SUBJECTS_DIR/${SUBNAME}/stats/rh.HCP-MMP1.stats -b \
			-a $SUBJECTS_DIR/${SUBNAME}/label/rh.HCP-MMP1.annot ${SUBNAME} rh white
	fi
done

#####
### Aggregate Regional Features Into Master Spreadsheet
#####

for SUBNAME in `ls $SUBJECTS_DIR | grep sub | grep -v ERROR | grep -v base` ; do
	if [[ -f `echo $SUBJECTS_DIR/$SUBNAME/stats/lh.HCP-MMP1.stats` ]] ; then 
		for HEMI in `echo lh rh` ; do
			cat $SUBJECTS_DIR/$SUBNAME/stats/${HEMI}.HCP-MMP1.stats | tail -n 180 > ${SUBNAME}_${HEMI}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$3'} | tr ' ' ',' | sed s@'ROI'@'ROI_area'@g >> ${SUBNAME}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$4'} | tr ' ' ',' | sed s@'ROI'@'ROI_volume'@g >> ${SUBNAME}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$5'} | tr ' ' ',' | sed s@'ROI'@'ROI_thickness'@g >> ${SUBNAME}.csv
		done
		cat ${SUBNAME}.csv | sed s@"^R_"@"rh_R_"@g | sed s@"^L_"@"lh_L_"@g > ${SUBNAME}_TEMP.csv ; mv ${SUBNAME}_TEMP.csv ${SUBNAME}.csv
		cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep -v ^#  | awk {'print $5,$4'} | tr ' ' ',' | sed s@'-'@'.'@g >> ${SUBNAME}.csv
		cat ${SUBNAME}.csv | grep -v hypointensities  | grep -v vessel | grep -v Optic | grep -v VentralDC | grep -v choroid > ${SUBNAME}_TEMP.csv
		cat ${SUBNAME}_TEMP.csv | sed s@'3rd'@'X3rd'@g | sed s@'4th'@'X4th'@g | grep -v '5th.Ventricle' | grep -v 'CSF' > ${SUBNAME}_TEMP1.csv
		echo lhCortexVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep lhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo rhCortexVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep rhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo lhCorticalWhiteMatterVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep lhCorticalWhiteMatterVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo rhCorticalWhiteMatterVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep rhCorticalWhiteMatterVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo SubCortGrayVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep SubCortGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo TotalGrayVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep TotalGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo SupraTentorialVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep SupraTentorialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo EstimatedTotalIntraCranialVol,`cat $SUBJECTS_DIR/$SUBNAME/stats/aseg.stats | grep EstimatedTotalIntraCranialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		cat ${SUBNAME}_TEMP1.csv | sed s@'-'@'.'@g > ${SUBNAME}_TEMP2.csv ; mkdir -p $PROJECT_DIR/freesurf_extract/Data 
		mv ${SUBNAME}_TEMP2.csv $PROJECT_DIR/freesurf_extract/Data/${SUBNAME}.csv ; rm ${SUBNAME}_*.csv
	fi
done

for SUBFILE in `ls $PROJECT_DIR/freesurf_extract/Data/sub-*.csv` ; do
	if [[ `cat $SUBFILE | wc -l` != 1118 ]] ; then
		echo "Incomplete Dataset Removed: `basename $SUBFILE` - `cat $SUBFILE | wc -l`" ; rm $SUBFILE
	fi
done

#####
### Download GTB Model From GitHub
#####

for SUBDIR in `echo $SUBJECTS_DIR/* | tr ' ' '\n' | grep sub | grep -v ERROR` ; do
	if [[ ! -f ${SUBDIR}/stats/avg.orig.nofix.euler ]] ; then
		RIGHT=`script -c "mris_euler_number ${SUBDIR}/surf/rh.orig.nofix" | grep ">" | cut -d ">" -f 1 | cut -d "=" -f 4 | cut -d " " -f 2`
		LEFT=`script -c "mris_euler_number ${SUBDIR}/surf/lh.orig.nofix" | grep ">" | cut -d ">" -f 1 | cut -d "=" -f 4 | cut -d " " -f 2`
		SUM=`echo $(( ${LEFT} + ${RIGHT} ))` ; AVG=`echo $(($SUM/2))` ; echo $AVG > ${SUBDIR}/stats/avg.orig.nofix.euler 
		rm typescript ; echo "Finishing Extraction for: `basename $SUBDIR`"
	fi
done

#####
### Download GTB Model From GitHub
#####

git clone https://github.com/RobertJirsaraie/brainage.git
mv brainage $SCRIPTS_DIR/GradientTreeBoost

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
