#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --job-name=DwiExtract
#############################

SUB=$1 
ROOT=/scratch/jirsaraie
HCD=/ceph/intradb/archive/CCF_HCD_STG/arc001
ATLASES=$ROOT/RobertJirsaraie/toolbox/bids_apps/atlases
DEPENDENCIES=$ROOT/RobertJirsaraie/toolbox/bids_apps/dependencies
PREPROCESSED=$HCD/$SUB/RESOURCES/Diffusion_preproc/$SUB/Diffusion/data

#####
### Install the MRTrix & AMICO Software Packages
#####

BIDDIR=$ROOT/study-HCD/bids/$SUB ; mkdir -p $BIDDIR/dwi
DWIDIR=$ROOT/study-HCD/apps/dti-fsl/$SUB ; mkdir -p $DWIDIR
NODDIDIR=$ROOT/study-HCD/apps/dti-noddi/$SUB ; mkdir -p $NODDIDIR
if [[ ! -f $DEPENDENCIES/mrtrix3.simg ]] ; then 
	singularity build $DEPENDENCIES/mrtrix3.simg docker://mrtrix3/mrtrix3
	pip install dmri-amico -U
fi

#####
### Copy Processed Data and Normalize Into Standardized Space
#####

cp $PREPROCESSED/nodif_brain_mask.nii.gz -sf $BIDDIR/dwi
cp $PREPROCESSED/nodif_brain.nii.gz -sf $BIDDIR/dwi
cp $PREPROCESSED/data.nii.gz -sf $BIDDIR/dwi
cp $PREPROCESSED/bvals -sf $BIDDIR/dwi
cp $PREPROCESSED/bvecs -sf $BIDDIR/dwi

cat $HCD/$SUB/RESOURCES/Diffusion_preproc/$SUB/Diffusion/QC/qc.json \
	| grep qc_mot_abs \
	| awk {'print $2'} \
	| sed s@','@''@g > $BIDDIR/dwi/qc_mot_abs.1D

cat $HCD/$SUB/RESOURCES/Diffusion_preproc/$SUB/Diffusion/QC/qc.json \
	| grep qc_mot_rel \
	| awk {'print $2'} \
	| sed s@','@''@g > $BIDDIR/dwi/qc_mot_rel.1D

#####
### Prepare the JHU Atlas 
#####

antsRegistrationSyNQuick.sh \
	-d 3 \
	-f $BIDDIR/dwi/nodif_brain.nii.gz \
	-m $BIDDIR/anat/T1w_restore_brain.nii.gz  \
	-o $BIDDIR/dwi/warp_T1w
antsApplyTransforms \
	-d 3 \
	-i $BIDDIR/anat/T1w_restore_brain.nii.gz \
	-r $BIDDIR/dwi/nodif_brain.nii.gz \
	-o $BIDDIR/dwi/anat_T1w.nii.gz \
	-t $BIDDIR/dwi/warp_T1w1Warp.nii.gz \
	-t $BIDDIR/dwi/warp_T1w0GenericAffine.mat \
	-n MultiLabel \
	-v 1
fast $BIDDIR/dwi/anat_T1w.nii.gz
rm $BIDDIR/dwi/warp_T1w* $BIDDIR/dwi/anat_T1w_*type.nii.gz $BIDDIR/dwi/*seg.nii.gz
antsRegistrationSyNQuick.sh \
	-d 3 \
	-f $BIDDIR/dwi/anat_T1w_pve_2.nii.gz \
	-m $ATLASES/atl-JHU_res-01_DWI.nii.gz \
	-o $BIDDIR/dwi/warp_atl
antsApplyTransforms \
	-d 3 \
	-i $ATLASES/atl-JHU_res-01_DWI.nii.gz \
	-r $BIDDIR/dwi/anat_T1w_pve_2.nii.gz \
	-o $BIDDIR/dwi/atl-JHU_DWI.nii.gz \
	-t $BIDDIR/dwi/warp_atl1Warp.nii.gz \
	-t $BIDDIR/dwi/warp_atl0GenericAffine.mat \
	-n MultiLabel \
	-v 1
rm $BIDDIR/dwi/warp_atl*

#####
### Compute DTI Maps using FSL
#####

dtifit \
	-k $BIDDIR/dwi/data.nii.gz \
	-m $BIDDIR/dwi/nodif_brain_mask.nii.gz  \
	-b $BIDDIR/dwi/bvals \
	-r $BIDDIR/dwi/bvecs \
	-o $DWIDIR/FIT 
rm `find $DWIDIR -type f | grep -v FA.nii | grep -v MD.nii`

#####
### Compute NODDI Maps using AMICO 
#####

python3 $ROOT/RobertJirsaraie/proj22-NeuroMap/Apps_DwiFit.py $SUB
mv `find $NODDIDIR | grep .nii.gz` $NODDIDIR ; rm -rf $NODDIDIR/$SUB
rm `find $NODDIDIR -type f | grep -v OD.nii | grep -v ICVF.nii`

#####
### Extract Regional Features of White Matter Microstructure
#####

for INDEX in `seq 1 48` ; do
	fslmaths $BIDDIR/dwi/atl-JHU_DWI.nii.gz -thr $INDEX -uthr $INDEX -bin $BIDDIR/dwi/TEMP_${INDEX}.nii.gz
	for MAP in `find $NODDIDIR $DWIDIR -type f | grep FIT_ | grep .nii.gz` ; do
		DIR=`dirname $MAP` ; LABEL=`basename $MAP | cut -d '.' -f1 | cut -d '_' -f2`
		fslmaths $MAP -mul $BIDDIR/dwi/TEMP_${INDEX}.nii.gz $DIR/${LABEL}_JHU-${INDEX}.nii.gz
		fslstats $DIR/${LABEL}_JHU-${INDEX}.nii.gz -M > $DIR/FEAT-${LABEL}_JHU-${INDEX}.1D
		rm $DIR/${LABEL}_JHU-${INDEX}.nii.gz
	done
	rm $BIDDIR/dwi/TEMP_${INDEX}.nii.gz
done
find $NODDIDIR -empty -type f -delete
find $DWIDIR -empty -type f -delete

#####
### Perform Deterministic Tractography using CAMINO
#####

<<FAILED_ANALYSIS

export PATH=${PATH}:${DEPENDENCIES}/camino/bin
export CAMINO_HEAP_SIZE=10000

fsl2scheme \
	-bvecfile $BIDDIR/dwi/bvecs \
	-bvalfile $BIDDIR/dwi/bvals \
	> $TRACDIR/bvals.scheme

image2voxel \
	-4dimage $BIDDIR/dwi/data.nii.gz \
	-outputfile $TRACDIR/dtfit.Bfloat 

dtfit \
	$TRACDIR/dtfit.Bfloat \
	$TRACDIR/bvals.scheme \
	-bgmask $BIDDIR/dwi/nodif_brain_mask.nii.gz \
	-outputfile $TRACDIR/dtfit.Bdouble 

track \
	-inputmodel dt \
	-seedfile $BIDDIR/dwi/warped_ATL.nii.gz \
	-inputfile $TRACDIR/dtfit.Bdouble \
	-tracker euler \
	-interpolator linear \
	-iterations 20 \
	-curvethresh 60 | \
	procstreamlines \
		-waypointfile $BIDDIR/dwi/warped_T1w_pve_2_dilate4.nii.gz \
		-exclusionfile $BIDDIR/dwi/warped_T1w_pve_0.nii.gz \
		-truncateinexclusion \
		-endpointfile $BIDDIR/dwi/warped_ATL.nii.gz \
		-outputfile $TRACDIR/tracts.Bdouble

#####
### Extract Microstructure Measures From Each White Matter Tract
#####

LABELS=`cat $ATLASES/atl-Gordon_Labels_Parcels.1D`
INDEX=`cat $ATLASES/atl-Gordon_Labels_Parcels.1D | sed s@'x'@' '@g | awk {'print $1'}`
paste <(printf %s "$INDEX") <(printf %s "$LABELS") > $TRACDIR/atl-Gordon_Parcels.tsv

for MAP in `find $ROOT/study-HCD/apps/dti-*/$SUB/FIT*.nii.gz` ; do
	OUTLABEL=$(echo `basename $MAP | cut -d '.' -f1`_)
	conmat \
		-inputfile $TRACDIR/tracts.Bdouble \
		-targetfile $BIDDIR/dwi/warped_ATL.nii.gz \
		-scalarfile $MAP \
		-tractstat mean \
		-targetnamefile $TRACDIR/atl-Gordon_Parcels.tsv \
		-outputroot $OUTLABEL
done

counttracts -inputfile $TRACDIR/tracts.Bdouble > $TRACDIR/tractCount.1D

#####
### Register the Anat & Atlas to DWI Space
#####

antsRegistrationSyNQuick.sh \
	-d 3 \
	-f $BIDDIR/dwi/nodif_brain.nii.gz \
	-m $BIDDIR/anat/T1w_restore_brain.nii.gz  \
	-o $BIDDIR/dwi/warp_T1w
antsApplyTransforms \
	-d 3 \
	-i $BIDDIR/anat/T1w_restore_brain.nii.gz \
	-r $BIDDIR/dwi/nodif_brain.nii.gz \
	-o $BIDDIR/dwi/warped_T1w.nii.gz \
	-t $BIDDIR/dwi/warp_T1w1Warp.nii.gz \
	-t $BIDDIR/dwi/warp_T1w0GenericAffine.mat \
	-n MultiLabel \
	-v 1
antsRegistrationSyNQuick.sh \
	-d 3 \
	-f $BIDDIR/dwi/nodif_brain.nii.gz \
	-m $ATLASES/atl-Gordon_T1w_Parcels.nii.gz  \
	-o $BIDDIR/dwi/warp_ATL
antsApplyTransforms \
	-d 3 \
	-i $ATLASES/atl-Gordon_T1w_Parcels.nii.gz  \
	-r $BIDDIR/dwi/nodif_brain.nii.gz \
	-o $BIDDIR/dwi/warped_ATL.nii.gz \
	-t $BIDDIR/dwi/warp_ATL1Warp.nii.gz \
	-t $BIDDIR/dwi/warp_ATL0GenericAffine.mat \
	-n MultiLabel \
	-v 1
rm -rf $BIDDIR/dwi/warp_*
fast $BIDDIR/dwi/warped_T1w.nii.gz
ImageMath 3 \
	$BIDDIR/dwi/warped_T1w_pve_2_dilate4.nii.gz GD \
	$BIDDIR/dwi/warped_T1w_pve_2.nii.gz 4

FAILED_ANALYSIS

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######