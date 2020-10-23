#!/usr/bin/env bash
#PBS -o /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stdout
#PBS -e /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stderr
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=12:00:00
#PBS -N MASS_PIPELINE
#PBS -q dque
#PBS -j oe 
source ~/ConfigEnv.sh
#####################

STUDY_DIR=/scratch/rjirsara/study-PDS
PROJECT_DIR=/scratch/rjirsara/projects/BrainAgeEval/DeepBrainNet
SCRIPTS_DIR=/scratch/rjirsara/RobertJirsaraie/projects/BrainAgeEval

git clone https://github.com/RobertJirsaraie/DeepBrainNet.git

python ${SCRIPTS_DIR}/DeepBrainNet/Script/Slicer.py ${SCRIPTS_DIR}/DeepBrainNet/Data/ ${SCRIPTS_DIR}/

python ${SCRIPTS_DIR}/DeepBrainNet/Script/Model_Test.py ${SCRIPTS_DIR}/Test ${SCRIPTS_DIR}/output.csv  ${SCRIPTS_DIR}/DeepBrainNet/Models/

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
