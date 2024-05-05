#!/usr/bin/env bash
#PBS -N AnatExtract
#PBS -l nodes=1:ppn=8,mem=100gb,walltime=24:00:00
#PBS -o /scratch/jirsaraie/qsubs/${PBS_JOBNAME}.log
#PBS -o /scratch/jirsaraie/qsubs/${PBS_JOBNAME}.error
#PBS -q dque
#PBS -j oe 
source ~/Anaconda_Env.sh
########################

ROOT=/scratch/jirsaraie
TOOLS=$ROOT/RobertJirsaraie/toolbox
export SUBJECTS_DIR=$ROOT/study-HCD/apps/freesurfer

#####
### Extract Regional Features from Each Atlas File
#####

for SUBNAME in `ls $ROOT/study-HCD/bids` ; do
	for ANNOT in `ls $TOOLS/bids_apps/atlases/atl-*.annot | grep 'atl-Gordon'` ; do
		ATLAS=`basename $ANNOT` 
		STATS=`echo $ATLAS | sed s@'.annot'@'.stats'@g`
		HEMI=`echo $ATLAS | cut -d '_' -f2 | cut -d '.' -f1`
		if [[ ! -f `echo $SUBJECTS_DIR/$SUBNAME/stats/${STATS}` ]] ; then
			mri_surf2surf \
				--srcsubject fsaverage \
				--trgsubject $SUBNAME \
				--hemi $HEMI \
				--sval-annot ${ANNOT} \
				--tval $SUBJECTS_DIR/$SUBNAME/label/${ATLAS}
			mris_anatomical_stats \
				-mgz -cortex $SUBJECTS_DIR/$SUBNAME/label/${HEMI}.cortex.label \
				-f $SUBJECTS_DIR/$SUBNAME/stats/${STATS} \
				-b -a $SUBJECTS_DIR/$SUBNAME/label/${ATLAS} $SUBNAME ${HEMI} white
		fi
	done
done

#####
### Aggregate Anatomical Features from the Glasser 2016 Atlas 
#####

<<FAILED_PARCELLATION
for SUBDIR in `echo $ROOT/study-HCD/apps/freesurfer/*` ; do
	if [[ ! -f `echo $SUBDIR/stats/lh.HCP-MMP1.stats` || ! -f `echo $SUBDIR/stats/glasser.csv` ]] ; then 
		SUBNAME=`basename $SUBDIR`
		for HEMI in `echo lh rh` ; do
			cat $SUBDIR/stats/${HEMI}.HCP-MMP1.stats | tail -n 180 > ${SUBNAME}_${HEMI}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$3'} | tr ' ' ',' | sed s@'ROI'@'ROI_area'@g >> ${SUBNAME}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$4'} | tr ' ' ',' | sed s@'ROI'@'ROI_volume'@g >> ${SUBNAME}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$5'} | tr ' ' ',' | sed s@'ROI'@'ROI_thickness'@g >> ${SUBNAME}.csv
			#cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$7'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_meancurv_parc"@g >> ${SUBNAME}.csv
			#cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$8'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_gauscurv_parc"@g >> ${SUBNAME}.csv
			#cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$9'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_foldindex_parc"@g >> ${SUBNAME}.csv
			#cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$10'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_curvindex_parc"@g >> ${SUBNAME}.csv
		done
		#cat $SUBDIR/stats/aseg.stats | grep -v ^#  | awk {'print $5,$4'} | tr ' ' ',' | sed s@'-'@'.'@g >> ${SUBNAME}.csv
		cat ${SUBNAME}.csv | sed s@"^R_"@"rh_R_"@g | sed s@"^L_"@"lh_L_"@g > ${SUBNAME}_TEMP.csv ; mv ${SUBNAME}_TEMP.csv ${SUBNAME}.csv
		cat $SUBDIR/stats/aseg.stats | grep -v ^#  | awk {'print $5,$4'} | tr ' ' ',' | sed s@'-'@'.'@g >> ${SUBNAME}.csv
		cat ${SUBNAME}.csv | grep -v hypointensities  | grep -v vessel | grep -v Optic | grep -v VentralDC | grep -v choroid > ${SUBNAME}_TEMP.csv
		cat ${SUBNAME}_TEMP.csv | sed s@'3rd'@'X3rd'@g | sed s@'4th'@'X4th'@g | grep -v '5th.Ventricle' | grep -v 'CSF' > ${SUBNAME}_TEMP1.csv
		echo lhCortexVol,`cat $SUBDIR/stats/aseg.stats | grep lhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo rhCortexVol,`cat $SUBDIR/stats/aseg.stats | grep rhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo lhCerebralWhiteMatter,`cat $SUBDIR/stats/aseg.stats | grep lhCerebralWhiteMatter | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo rhCerebralWhiteMatter,`cat $SUBDIR/stats/aseg.stats | grep rhCerebralWhiteMatter | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo SubCortGrayVol,`cat $SUBDIR/stats/aseg.stats | grep SubCortGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo TotalGrayVol,`cat $SUBDIR/stats/aseg.stats | grep TotalGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo SupraTentorialVol,`cat $SUBDIR/stats/aseg.stats | grep SupraTentorialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo EstimatedTotalIntraCranialVol,`cat $SUBDIR/stats/aseg.stats | grep EstimatedTotalIntraCranialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		cat ${SUBNAME}_TEMP1.csv | sed s@'-'@'.'@g > ${SUBNAME}_TEMP2.csv
		mv ${SUBNAME}_TEMP2.csv $SUBDIR/stats/glasser.csv 
		rm ${SUBNAME}*.csv
	fi
done

for SUBFILE in `ls $ROOT/study-HCD/apps/freesurfer/*/stats/glasser.csv` ; do
	if [[ `cat $SUBFILE | wc -l` != 1118 ]] ; then
		echo "Incomplete Dataset Removed: `basename $SUBFILE` - `cat $SUBFILE | wc -l`" ; rm $SUBFILE
	fi
done
FAILED_PARCELLATION

#####
### Aggregate Anatomical Features from the Gordon Atlas 
#####

for SUBDIR in `echo $ROOT/study-HCD/apps/freesurfer/HCD*` ; do
if [[ -f `echo $SUBDIR/stats/atl-Gordon_rh.HCP-MMP1.stats` && ! -f `echo $SUBDIR/stats/gordon.csv` ]] ; then 
		SUBNAME=`basename $SUBDIR`
		for HEMI in `echo lh rh` ; do
			if [ $HEMI == 'lh' ] ; then
				cat $SUBDIR/stats/atl-Gordon_${HEMI}.HCP-MMP1.stats | tail -n 161 > ${SUBNAME}_${HEMI}.csv
			else 
				cat $SUBDIR/stats/atl-Gordon_${HEMI}.HCP-MMP1.stats | tail -n 172 > ${SUBNAME}_${HEMI}.csv
			fi
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$3'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_area_parcel"@g >> ${SUBNAME}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$4'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_volume_parcel"@g >> ${SUBNAME}.csv
			cat ${SUBNAME}_${HEMI}.csv | awk {'print $1,$5'} | tr ' ' ',' | sed s@'Parcel'@"${HEMI}_thickness_parcel"@g >> ${SUBNAME}.csv
		done
		cat $SUBDIR/stats/aseg.stats | grep -v ^#  | awk {'print $5,$4'} | tr ' ' ',' | sed s@'-'@'.'@g >> ${SUBNAME}.csv
		cat ${SUBNAME}.csv | grep -v hypointensities  | grep -v vessel | grep -v Optic | grep -v VentralDC | grep -v choroid > ${SUBNAME}_TEMP.csv
		cat ${SUBNAME}_TEMP.csv | sed s@'3rd'@'X3rd'@g | sed s@'4th'@'X4th'@g | grep -v '5th.Ventricle' | grep -v 'CSF' > ${SUBNAME}_TEMP1.csv
		echo lhCortexVol,`cat $SUBDIR/stats/aseg.stats | grep lhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo rhCortexVol,`cat $SUBDIR/stats/aseg.stats | grep rhCortexVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo lhCerebralWhiteMatter,`cat $SUBDIR/stats/aseg.stats | grep lhCerebralWhiteMatter | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo rhCerebralWhiteMatter,`cat $SUBDIR/stats/aseg.stats | grep rhCerebralWhiteMatter | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo SubCortGrayVol,`cat $SUBDIR/stats/aseg.stats | grep SubCortGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo TotalGrayVol,`cat $SUBDIR/stats/aseg.stats | grep TotalGrayVol | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo SupraTentorialVol,`cat $SUBDIR/stats/aseg.stats | grep SupraTentorialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		echo EstimatedTotalIntraCranialVol,`cat $SUBDIR/stats/aseg.stats | grep EstimatedTotalIntraCranialVol, | cut -d ',' -f4 | sed s@' '@''@g` >> ${SUBNAME}_TEMP1.csv
		cat ${SUBNAME}_TEMP1.csv | sed s@'-'@'.'@g > ${SUBNAME}_TEMP2.csv
		mv ${SUBNAME}_TEMP2.csv $SUBDIR/stats/gordon.csv 
		rm ${SUBNAME}*.csv
	fi
done

for SUBFILE in `ls $ROOT/study-HCD/apps/freesurfer/*/stats/gordon.csv` ; do
	COUNT=`cat $SUBFILE | wc -l`
	if [ $COUNT != 2369 ] ; then
		if [ $COUNT != 1037 ] ; then
			echo "Incomplete Dataset Removed: `echo $SUBFILE` - `cat $SUBFILE | wc -l`" ; rm $SUBFILE
		fi
	fi
done 

#####
### Extract the Euler Number For Each Subject
#####

for SUBDIR in `echo $ROOT/study-HCD/apps/freesurfer/* | tr ' ' '\n' | grep -v fsaverage` ; do
	if [[ ! -f ${SUBDIR}/stats/avg.orig.nofix.euler ]] ; then
		RIGHT=`script -c "mris_euler_number ${SUBDIR}/surf/rh.orig.nofix" | grep ">" | cut -d ">" -f 1 | cut -d "=" -f 4 | cut -d " " -f 2`
		LEFT=`script -c "mris_euler_number ${SUBDIR}/surf/lh.orig.nofix" | grep ">" | cut -d ">" -f 1 | cut -d "=" -f 4 | cut -d " " -f 2`
		SUM=`echo $(( ${LEFT} + ${RIGHT} ))` ; AVG=`echo $(($SUM/2))` ; echo $AVG > ${SUBDIR}/stats/avg.orig.nofix.euler 
		rm typescript ; echo "Finishing Extraction for: `basename $SUBDIR`"
	fi
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######