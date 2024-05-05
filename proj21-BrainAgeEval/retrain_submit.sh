#!/usr/bin/env bash
###################

DataFreeze=/scratch/jirsaraie/proj20-BrainAgeEval/study-HCP/Analysis/n789_DataFreeze_20210305.csv
ModelDir=/scratch/jirsaraie/RobertJirsaraie/proj20-BrainAgeEval/DeepBrainNet/Models
ImgDir=/scratch/jirsaraie/proj20-BrainAgeEval/study-HCP/DBNInput

######
### Create Cohort File
######

if [[ ! -f `echo ${ImgDir}/Retrain.csv` ]] ; then
	echo "sub,path" > ${ImgDir}/Retrain.csv
	for DIR in `echo ${ImgDir}/*NormCor_brain*.jpg` ; do
		SUB=`basename ${DIR}  | sed s@_Norm@#@g | cut -d '#' -f1`
		echo ${SUB},`basename ${DIR}` >> ${ImgDir}/Retrain.csv
	done
fi

######
### Create Basic Structure of the Submission Scripts
######

for Split in `seq 1 1 5` ; do
	
	ScriptDir=/scratch/jirsaraie/RobertJirsaraie/proj20-BrainAgeEval/DeepBrainNet/Script
	QUEUEFULL=${ScriptDir}/submit_retrain_f${Split}.sh ; QUEUEHALF=${ScriptDir}/submit_retrain_h${Split}.sh
	echo '#!/bin/bash' > ${QUEUEFULL}
	echo "#PBS -o ${ScriptDir}/submit_retrain_f${Split}.stdout" >> ${QUEUEFULL}
	echo "#PBS -e ${ScriptDir}/submit_retrain_f${Split}.stderr" >> ${QUEUEFULL}
	echo "#PBS -l nodes=1:ppn=16,mem=16gb,walltime=168:00:00" >> ${QUEUEFULL}
	echo "#PBS -N DBNtrainF${Split}" >> ${QUEUEFULL}
	echo "#PBS -q dque" >> ${QUEUEFULL}
	echo "#PBS -j oe" >> ${QUEUEFULL}
	echo "conda activate tensorflow" >> ${QUEUEFULL}
	echo "#########################" >> ${QUEUEFULL}
	echo "" >> ${QUEUEFULL}

	cp $QUEUEFULL $QUEUEHALF ; cat $QUEUEHALF | sed s@'retrain_f'@'_retrain_h'@g > temp.sh ; mv temp.sh $QUEUEHALF 

	echo 'python3 /scratch/jirsaraie/RobertJirsaraie/proj20-BrainAgeEval/retrain_model.py \'  >> ${QUEUEFULL}
	echo '"'${DataFreeze}'"' '\' >> ${QUEUEFULL}
	echo '"'${ImgDir}'"' '\' >>     ${QUEUEFULL}
	echo '"'${ModelDir}'"' '\' >>   ${QUEUEFULL}
	echo "${Split} '"'full'"' " >>  ${QUEUEFULL}
	echo "" >> ${QUEUEFULL} ; chmod ug+wrx ${QUEUEFULL}
	qsub ${QUEUEFULL}

	echo 'python3 /scratch/jirsaraie/RobertJirsaraie/proj20-BrainAgeEval/retrain_model.py \'  >> ${QUEUEHALF}
	echo '"'${DataFreeze}'"' '\' >> ${QUEUEHALF}
	echo '"'${ImgDir}'"' '\' >>     ${QUEUEHALF}
	echo '"'${ModelDir}'"' '\' >>   ${QUEUEHALF}
	echo "${Split} '"'half'"' " >>  ${QUEUEHALF}
	echo "" >> ${QUEUEHALF} ; chmod ug+wrx ${QUEUEHALF}
	qsub ${QUEUEHALF}

done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
