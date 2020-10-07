#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=4:00:00
#PBS -q dque
source ~/ConfigEnv.sh
#####################

DIR_PROJECT=$1
ZIP_DOWNLOAD=$2

wget https://cnda.wustl.edu/xapi/archive/download/$ZIP_DOWNLOAD/zip --user=${USER} --password=${PASSWORD_NIL} --no-check-certificate -O $DIR_PROJECT/${ZIP_DOWNLOAD}.zip

unzip -u $DIR_PROJECT/${ZIP_DOWNLOAD}.zip -d $DIR_PROJECT

rm $DIR_PROJECT/${ZIP_DOWNLOAD}.zip ; rm XnatDL.e* XnatDL.o* ; chmod -R 750 $DIR_PROJECT

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
