#!/bin/bash
#SBATCH --job-name="multhelical"
#SBATCH --export=ALL
#SBATCH --exclusive
## number of nodes
#SBATCH -N 9
## number of cores
#SBATCH -n 216

cd $SLURM_SUBMIT_DIR

# These commands / variables should be modified to load the appropriate environment for the user's cluster.
# EMAN1, EMAN2/SPARX, Open MPI, and hsearch_lorentz (from Egelman) are required.
# The paths to eman1 and hsearch_lorentz are hardcoded in hfunctions.py :(
# refine.py is the script included in the protocol
HELICAL_REFINE_HOME=/rugpfs/fs0/home/hpcsoft/soft/helical_refine_0.0.1
source $HELICAL_REFINE_HOME/eman2_2.1/eman2.bashrc
export LD_LIBRARY_PATH=$HELICAL_REFINE_HOME/openmpi-3.0.2/lib:$LD_LIBRARY_PATH
export PATH=$HELICAL_REFINE_HOME/openmpi-3.0.2/bin:$PATH

PATH_TO_REFINE_SCRIPT="$HELICAL_REFINE_HOME/helical_refine/refine.py"

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

mpirun -np 216 $PATH_TO_REFINE_SCRIPT start.hdf init.hdf refine_dir --ou=120 --rs=1 --xr='16 8 4' --ts='4 2 1' --delta='4 3 2' --an='-1 30 10' --snr=0.08 --maxit=2 --ref_a=S --cutoff=1.5 --MPI --olmask=70 --ilmask=75 --lmask=280 --protos='11 12 13 14 15' --hsearch='75.0 190.0' --hpars='-32.503 12.289 -29.893 10.254 -27.678 9.495 -25.604 9.351 -23.82 10.984' --oplane=10 --recon_pad=2 --findseam --full_output --sort > refine.log

