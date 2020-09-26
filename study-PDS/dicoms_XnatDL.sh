#!/bin/sh
#PBS -N XnatDL
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=2gb
#PBS -q dque
############

DOWNLOAD=$1

wget https://cnda.wustl.edu/xapi/archive/download/${DOWNLOAD}/zip --user=rjirsara --password=Rjj219619! --no-check-certificate -O /scratch/rjirsara/study-PDS/dicoms/${DOWNLOAD}.zip

echo /scratch/rjirsara/study-PDS/dicoms/${DOWNLOAD}.zip
unzip -u /scratch/rjirsara/study-PDS/dicoms/${DOWNLOAD} -d /scratch/rjirsara/study-PDS/dicoms/

#####⚡⚡⚡⚡⚡⚡##############################⚡⚡⚡⚡⚡⚡#############################⚡⚡⚡⚡⚡⚡####
###         ⚡   ⚡   ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡   ⚡   ⚡       ###
#####⚡⚡⚡⚡⚡⚡##############################⚡⚡⚡⚡⚡⚡#############################⚡⚡⚡⚡⚡⚡####
