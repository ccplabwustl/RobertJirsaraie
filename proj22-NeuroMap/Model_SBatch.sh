#!/bin/bash
###########

TODAY=$(date +%Y%m%d)
PROJ=/scratch/jirsaraie/proj22-NeuroMap
STORE=$PROJ/data/featmap_2319x5_20230531.csv
FEATURES=`cat $STORE | tr ',' ' ' | awk {'print $2'} | tail -n +2 | sort | uniq`
NETWORKS=`cat $STORE | tr ',' ' ' | awk {'print $3'} | tail -n +2 | grep -v NA | sort | uniq`
for ITERATION in `echo FEAT-Multimod NET-Global $FEATURES $NETWORKS | sed s@'"'@''@g` ; do
	if [[ `echo $FEATURES | tr ' ' '\n' | grep $ITERATION | wc -l` == 1 ]] ; then
		ITERATION=`echo FEAT-${ITERATION}`
	fi
	if [[ `echo $NETWORKS | tr ' ' '\n' | grep $ITERATION | wc -l` == 1 ]] ; then
		ITERATION=`echo NET-${ITERATION}`
	fi
	if [ ! -d `echo $PROJ/analysis/${TODAY}_TEMP/${ITERATION}` ]; then
		echo '#!/bin/bash' > ${ITERATION}.sh
		echo "#SBATCH -J ${ITERATION}" >> ${ITERATION}.sh
		echo '#SBATCH -N 1' >> ${ITERATION}.sh
		echo '#SBATCH --gres gpu' >> ${ITERATION}.sh
		echo '#SBATCH -t 48:00:00' >> ${ITERATION}.sh 
		echo '###################' >> ${ITERATION}.sh
		echo '' >> ${ITERATION}.sh
		echo "python3 /scratch/jirsaraie/RobertJirsaraie/proj22-NeuroMap/Model_Optimize.py ${ITERATION} ${TODAY}" >> ${ITERATION}.sh 
		echo "Submitting Job For ${ITERATION}" ; sbatch ${ITERATION}.sh ; rm ${ITERATION}.sh
	fi
done

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######