#!/usr/bin/env python
#PBS -o /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stdout
#PBS -e /scratch/rjirsara.submitted_jobs/${PBS_JOBNAME}_${PBS_JOBID}.stderr
#PBS -l nodes=1:ppn=8,mem=16gb,walltime=1:00:00
#PBS -N COMBAT_GAM
#PBS -q dque
#PBS -j oe 
#source ~/ConfigEnv.sh
#####################
#pip install git+https://github.com/statsmodels/statsmodels
#pip install git+https://github.com/rpomponio/neuroHarmonize

from neuroHarmonize import harmonizationLearn, harmonizationApply, loadHarmonizationModel
import pandas as pd
import numpy as np

######
### Apply GAM COMBAT to the Cross-sectional Freesurfer Processed Data Stream
######

data = pd.read_csv('/scratch/rjirsara/proj20-BrainAgeEval/analysis/n829_FS-Cross_20201231.csv')
ids = data[['sub','ses']]
data.pop("sub") 
data.pop("ses")
data.pop("sex")
data = np.array(data)

covars = pd.read_csv('/scratch/rjirsara/proj20-BrainAgeEval/analysis/n834_DataFreeze_20201231.csv')
covars = pd.merge(ids, covars,  how='left', left_on=['sub','ses'], right_on = ['sub','ses'])
covars[["SITE"]]=covars[["ses"]].replace([1,2,3],1) 
covars[["SITE"]]=covars[["SITE"]].replace([4,5],2)
covars = covars[['age','sex','SITE']]
covars[['age','sex','SITE']].apply(pd.to_numeric, errors='coerce', axis=1)

model, data_adj = harmonizationLearn(data, covars,smooth_terms=['age'])
CROSS = pd.concat([ids,covars], axis=1) ; CROSS = pd.concat([CROSS,pd.DataFrame(data_adj)], axis=1)
pd.DataFrame(CROSS).to_csv('/scratch/rjirsara/proj20-BrainAgeEval/analysis/n829_COMBATxGAM-Cross_20210101.csv', index=False)

######
### Apply GAM COMBAT to the Longitudinal Freesurfer Processed Data Stream
######

data = pd.read_csv('/scratch/rjirsara/proj20-BrainAgeEval/analysis/n828_FS-Long_20201231.csv')
ids = data[['sub','ses']]
data.pop("sub") 
data.pop("ses")
data.pop("sex")
data = np.array(data)

covars = pd.read_csv('/scratch/rjirsara/proj20-BrainAgeEval/analysis/n834_DataFreeze_20201231.csv')
covars = pd.merge(ids, covars,  how='left', left_on=['sub','ses'], right_on = ['sub','ses'])
covars[["SITE"]]=covars[["ses"]].replace([1,2,3],1) 
covars[["SITE"]]=covars[["SITE"]].replace([4,5],2)
covars = covars[['age','sex','SITE']]
covars[['age','sex','SITE']].apply(pd.to_numeric, errors='coerce', axis=1)

model, data_adj = harmonizationLearn(data, covars,smooth_terms=['age'])
LONG = pd.concat([ids,covars], axis=1) ; LONG = pd.concat([LONG,pd.DataFrame(data_adj)], axis=1)
pd.DataFrame(LONG).to_csv('/scratch/rjirsara/proj20-BrainAgeEval/analysis/n828_COMBATxGAM-Long_20210101.csv', index=False)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
