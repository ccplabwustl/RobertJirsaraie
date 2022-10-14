#!/usr/bin/env bash
###################

OUTDIR=/scratch/jirsaraie/collaborations/19-YassaLab
Realease=/scratch/jirsaraie/study-HCD/data/qa_1645x13_20220822.csv

#####
### Reformat files to BIDs specifications including relevant metadata
#####

sub-0330122/ses-1

for SUBJECT in `cat $Realease | tr ',' ' ' | awk {'print $1,$13'} | grep "Include" | awk {'print $1'} | sed s@'"'@''@g` ; do

	SUBDIR=/ceph/intradb/archive/CinaB/CCF_HCD_STG/${SUBJECT}_MR
	SUB=`echo $SUBJECT | cut -d '_' -f1 | sed s@'HCD'@''@g`
	SES=`echo $SUBJECT | cut -d '_' -f2 | sed s@V@@g`

	if [ ! -d $OUTDIR/bids/sub-${SUB}/ses-${SES}/func ] ; then
		echo $SUBDIR

		#Anatomical
		mkdir -p $OUTDIR/bids/sub-${SUB}/ses-${SES}/anat

		cp $SUBDIR/unprocessed/T1w_MPR_vNav_4e_e1e2_mean/*_T1w_MPR_vNav_4e_e1e2_mean.nii.gz \
			$OUTDIR/bids/sub-${SUB}/ses-${SES}/anat/sub-${SUB}_ses-${SES}_T1w.nii.gz

		cp $SUBDIR/unprocessed/T1w_MPR_vNav_4e_e1e2_mean/*_T1w_MPR_vNav_4e_e1e2_mean.json \
			$OUTDIR/bids/sub-${SUB}/ses-${SES}/anat/sub-${SUB}_ses-${SES}_T1w.json

		cp $SUBDIR/unprocessed/T2w_SPC_vNav/*_T2w_SPC_vNav.nii.gz \
			$OUTDIR/bids/sub-${SUB}/ses-${SES}/anat/sub-${SUB}_ses-${SES}_T2w.nii.gz

		cp $SUBDIR/unprocessed/T2w_SPC_vNav/*_T2w_SPC_vNav.json \
			$OUTDIR/bids/sub-${SUB}/ses-${SES}/anat/sub-${SUB}_ses-${SES}_T2w.json

		#Functional
		mkdir -p $OUTDIR/bids/sub-${SUB}/ses-${SES}/func
		for RUN in `seq 1 1 $(echo $SUBDIR/unprocessed/rfMRI_REST* | wc -w)` ; do
			FUNCDIR=`echo $SUBDIR/unprocessed/rfMRI_REST* | tr ' ' '\n' | head -n${RUN} | tail -n1`
			DIR=`basename $FUNCDIR | cut -d '_' -f3`

			sed '16i\
				\ "TaskName": "rest",' $FUNCDIR/*_rfMRI_REST*${DIR}.json > \
				$OUTDIR/bids/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-rest_run-${RUN}_bold.json

			cp $FUNCDIR/*_rfMRI_REST*_${DIR}.nii.gz \
				$OUTDIR/bids/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-rest_run-${RUN}_bold.nii.gz

			cp $FUNCDIR/*_rfMRI_REST*_${DIR}_SBRef.nii.gz \
				$OUTDIR/bids/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-rest_run-${RUN}_sbref.nii.gz

			### Excluded because we do not need fmaps for the ConteCenter processing stream 

			#cp $FUNCDIR/*_SpinEchoFieldMap*_AP.nii.gz  \
			#	$OUTDIR/bids/sub-${SUB}/ses-${SES}/fmap/sub-${SUB}_ses-${SES}_task-REST_dir-AP_run-${RUN}_epi.nii.gz

			#cp $FUNCDIR/*_SpinEchoFieldMap*_AP.json \
			#	$OUTDIR/bids/sub-${SUB}/ses-${SES}/fmap/sub-${SUB}_ses-${SES}_task-REST_dir-AP_run-${RUN}_epi.json

			#cp $FUNCDIR/*_SpinEchoFieldMap*_PA.nii.gz \
			#	$OUTDIR/bids/sub-${SUB}/ses-${SES}/fmap/sub-${SUB}_ses-${SES}_task-REST_dir-PA_run-${RUN}_epi.nii.gz

			#cp $FUNCDIR/*_SpinEchoFieldMap*_PA.json \
			#	$OUTDIR/bids/sub-${SUB}/ses-${SES}/fmap/sub-${SUB}_ses-${SES}_task-REST_dir-PA_run-${RUN}_epi.json
		done
	fi

	if [ ! -d $OUTDIR/proc/sub-${SUB}/ses-${SES} ] ; then

		mkdir -p $OUTDIR/proc/sub-${SUB}/ses-${SES}
		PROCDIR=/scratch/jirsaraie/study-HCD/bids/HCD${SUB}_V${SES}_MR/func

		cp $SUBDIR/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL.nii.gz \
			$OUTDIR/proc/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-rest_preproc.nii.gz

		cp $PROCDIR/HCD${SUB}_V${SES}_MR_fully_denoised.nii.gz  \
			$OUTDIR/proc/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-rest_denoised.nii.gz

		cp $PROCDIR/FDfilt4.txt \
			$OUTDIR/proc/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-rest_preproc_FD.txt

		cp $PROCDIR/FDfilt4_censored.txt \
			$OUTDIR/proc/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-rest_denoised_FD.txt

	fi
done

bids-validator $OUTDIR/bids --ignore-warnings > $OUTDIR/bids_verification.txt

#####
### Upload the BIDs-Formatted Data into Flywheel
#####

fw login $API_KEY

fw import folder $OUTDIR/bids --group yassalab --project HCP-D3-BIDs --no-audit-log --skip-existing

fw import folder $OUTDIR/proc --group yassalab --project HCP-D3-Proc --no-audit-log --skip-existing

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####          ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡          ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######