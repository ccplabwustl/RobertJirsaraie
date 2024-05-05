#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################

import amico, os, sys

######
### Define Subject Paths
######

SUB=sys.argv[1]
BIDS=f"/scratch/jirsaraie/study-HCD/bids/{SUB}/dwi"
APPS=f"/scratch/jirsaraie/study-HCD/apps/dti-noddi/{SUB}"
os.chdir(APPS)

######
### Load Diffusion Data
######

amico.setup() 
BVAL=f"{BIDS}/bvals" 
BVEC=f"{BIDS}/bvecs"
amico.util.fsl2scheme(BVAL, BVEC) 
ae = amico.Evaluation(f"{APPS}/{SUB}", SUB)
ae.load_data(f"{BIDS}/data.nii.gz", f"{BIDS}/bvals.scheme", mask_filename=f"{BIDS}/nodif_brain_mask.nii.gz", b0_thr=5)

######
### Fit Diffusion Model
######

ae.set_model("NODDI")
ae.generate_kernels(regenerate=False)
ae.load_kernels()
ae.fit()

######
### Save Spatial Maps
######

ae.save_results(APPS)

########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######
####           ⚡     ⚡    ⚡   ⚡   ⚡   ⚡  ⚡  ⚡  ⚡    ⚡  ⚡  ⚡  ⚡   ⚡   ⚡   ⚡    ⚡     ⚡         ####
########⚡⚡⚡⚡⚡⚡#################################⚡⚡⚡⚡⚡⚡################################⚡⚡⚡⚡⚡⚡#######