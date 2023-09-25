#!/usr/bin/env bash
###################

ROOT=/scratch/jirsaraie
HCD=/ceph/intradb/archive/CinaB/CCF_HCD_STG
PREPROC=/ceph/intradb/archive/CCF_HCD_STG/arc001

#####
### Create SymLinks or Process Data to Project-Specific Directory
#####

for SUB in `ls $HCD/*/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_Runs.csv | cut -d '/' -f7` ; do

	###Anatomical PreProc
	mkdir -p $ROOT/study-HCD/bids/$SUB/anat $ROOT/.qsubs
	if [[ ! -f $ROOT/study-HCD/bids/$SUB/anat/T1w.nii.gz ]] ; then
		echo "${SUB} - COPY ANAT" 
		mkdir -p $ROOT/study-HCD/bids/$SUB/anat 
		cp -rsu $HCD/$SUB/MNINonLinear/T*.nii.gz $ROOT/study-HCD/bids/$SUB/anat
	fi

	###Anatomical FreesurferProc
	if [[ ! -d $ROOT/study-HCD/apps/freesurfer/$SUB ]] ; then
		echo "${SUB} - COPY FREESURFER" 
		mkdir -p $ROOT/study-HCD/apps/freesurfer
		cp -rsu $HCD/$SUB/T1w/$SUB $ROOT/study-HCD/apps/freesurfer/
	fi

	###Functional Denoising
	mkdir -p $ROOT/study-HCD/bids/$SUB/func $ROOT/.qsubs
	PREPROC=$HCD/$SUB/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL.nii.gz
	DENOISED=$ROOT/study-HCD/bids/$SUB/func/${SUB}_fully_denoised.nii.gz
	if [[ -f $PREPROC && ! -f $DENOISED ]] ; then
		echo "${SUB} - DENOISED fMRI" 
		cp -rsu $HCD/$SUB/MNINonLinear/brainmask_fs.2.nii.gz $ROOT/study-HCD/bids/$SUB/func/${SUB}_rfMRI_mask.nii.gz
		CALL='matlab -nodisplay -nosplash -nodesktop -r "run('\''$ROOT/.qsubs/${SUB}_denoise_bold'\'');exit;"'
		cat $ROOT/RobertJirsaraie/proj22-NeuroMap/Apps_FuncProc.m | \
		sed s@"SUBID"@"${SUB}"@g > $ROOT/.qsubs/${SUB}_denoise_bold.m
		echo '#!/usr/bin/env bash' > $ROOT/.qsubs/${SUB}_denoise_bold.sh
		echo 'source ~/Anaconda_Env.sh' >> $ROOT/.qsubs/${SUB}_denoise_bold.sh
		echo "ROOT=$ROOT ; SUB=$SUB" >> $ROOT/.qsubs/${SUB}_denoise_bold.sh
		echo $CALL >> $ROOT/.qsubs/${SUB}_denoise_bold.sh
		chmod u+wrx $ROOT/.qsubs/${SUB}_denoise_bold.sh
		qsub -N ${SUB}_PRE \
			-o $ROOT/.qsubs/${SUB}_denoise_bold.log \
			-e $ROOT/.qsubs/${SUB}_denoise_bold.error \
			-l nodes=1:ppn=8,mem=500gb,walltime=24:00:00 \
			-q dque \
			-j oe $ROOT/.qsubs/${SUB}_denoise_bold.sh
	fi
 	
	###Functional Feature Extraction
	if [[ -f $PREPROC && ! -f $ROOT/study-HCD/apps/xcp-dualreg/${SUB}/stats/node-1_cohesion.txt ]] ; then
		echo "${SUB} - EXTRACT fMRI" ; cd $ROOT/study-HCD/apps 
		mkdir -p xcp-alff xcp-dualreg xcp-fcon xcp-reho 
		cat $ROOT/RobertJirsaraie/proj22-NeuroMap/Apps_FuncExtract.sh | \
		sed s@"SUBID"@"${SUB}"@g > $ROOT/.qsubs/${SUB}_extract_bold.sh
		echo '#!/usr/bin/env bash' > $ROOT/.qsubs/${SUB}_temp_bold.sh
		echo 'source ~/Anaconda_Env.sh' >> $ROOT/.qsubs/${SUB}_temp_bold.sh
		echo $ROOT/.qsubs/${SUB}_extract_bold.sh >> $ROOT/.qsubs/${SUB}_temp_bold.sh
		chmod u+wrx $ROOT/.qsubs/${SUB}_extract_bold.sh
		qsub -N ${SUB}_EXTRACT \
			-o $ROOT/.qsubs/${SUB}_extract_bold.log \
			-e $ROOT/.qsubs/${SUB}_extract_bold.error \
			-l nodes=1:ppn=8,mem=100gb,walltime=24:00:00 \
			-q dque \
			-j oe $ROOT/.qsubs/${SUB}_temp_bold.sh
		rm $ROOT/.qsubs/${SUB}_temp_bold.sh
	fi

	###Diffusion Feature Extraction
	NODDI=$ROOT/study-HCD/apps/dti-noddi/${SUB}/FIT_ICVF.nii.gz
	RAW=$PREPROC/$SUB/RESOURCES/Diffusion_preproc/$SUB/Diffusion/data/data.nii.gz
	if [[ -f $RAW && ! -f $NODDI ]] ; then
		echo "${SUB} - EXTRACT DWI"
		LOG=$ROOT/.qsubs/${SUB}_extract_dwi.log ; mkdir -p $ROOT/.qsubs
		sbatch $ROOT/RobertJirsaraie/proj22-NeuroMap/Apps_DwiExtract.sh ${SUB} &> $LOG 
	fi
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####          ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡          ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######