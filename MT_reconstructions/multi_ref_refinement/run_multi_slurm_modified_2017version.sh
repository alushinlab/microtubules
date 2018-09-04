#!/bin/bash
#SBATCH --job-name="multhelical"
#SBATCH --export=ALL

cd $SLURM_SUBMIT_DIR

# These commands / variables should be modified to load the appropriate environment for the user's cluster.
# EMAN1, EMAN2/SPARX, Open MPI, and hsearch_lorentz (from Egelman) are required.
# refine.py is the script included in the protocol

PATH_TO_REFINE_SCRIPT="/ru-auth/local/home/sti/for_jove/scripts/EMAN2/EMAN2_actin2/refine.py"

# Definition of paramters:
# np : Number of mpi processes to launch
# ou : Outer radius of the reconstruction in X
# olmask : radius of outer mask to apply along helix in pixels (default = 15)
# lmask : length of mask to apply along helix in Angstroms (default=280)
# protos : Number of protofilaments, relevant for microtubule reconstructions.  For single-start helices, set to 1
# hpars : Initial guess of helical twist (degrees) and rise (Angstroms) for each model.  Left-handed twist is negative.
# hsearch : Inner and outer radius for helical search by hsearch_lorentz, in Angstroms
# 
# optional parameter:
# ilmask : radius of inner mask to apply along helix in Angstroms (default=0, change for hollow tubes like TMV or microtubules)

#execute

mpirun.eman2 -np 216 $PATH_TO_REFINE_SCRIPT start.hdf init.hdf refine --ou=120 --rs=1 --xr='16 8 4' --ts='4 2 1' --delta='4 3 2' --an='-1 30 10' --snr=0.08 --maxit=2 --ref_a=S --cutoff=1.5 --MPI --olmask=70 --ilmask=75 --lmask=280 --protos='13 14 15 15' --hsearch='75.0 190.0' --hpars='-27.617 10.194 -25.713 8.751 -24.009 8.74 -23.811 10.859' --oplane=10 --recon_pad=2 --findseam --full_output --sort > refine.log
