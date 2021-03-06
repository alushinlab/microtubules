#!/bin/bash
#PBS -N helical

cd $PBS_O_WORKDIR

# These commands / variables should be modified to load the appropriate environment for the user's cluster.
# EMAN1, EMAN2/SPARX, Open MPI, and hsearch_lorentz (from Egelman) are required.
# refine.py is the script included in 

MPI_COMMAND = "/usr/local/apps/EMAN2/openmpi/1.6.5/gnu/eth/bin/mpirun" 
PATH_TO_REFINE_SCRIPT = "/data/Alushinlab/lab_scripts/EMAN2/EMAN2_actin2/refine.py"

source /data/Alushinlab/environment/.bashrc
source /usr/local/EMAN1-1.9/eman.bashrc
module load eman2

# Definition of paramters:
# np : Number of mpi processes to launch
# ou : Outer radius of the reconstruction in X
# olmask : radius of outer mask to apply along helix in pixels (default = 15)
# lmask : length of mask to apply along helix in Angstroms (default=280)
# protos : Number of protofilaments, relevant for microtubule reconstructions.  For single-start helices, set to 1
# hpars : Initial guess of helical twist (degrees) and rise (Angstroms).  Left-handed twist is negative
# hsearch : Inner and outer radius for helical search by hsearch_lorentz, in Angstroms
# 
# optional parameter:
# ilmask : radius of inner mask to apply along helix in Angstroms (default=0, change for hollow tubes like TMV or microtubules)

#execute

$MPI_COMMAND -machinefile $PBS_NODEFILE -np 32 $PATH_TO_REFINE_SCRIPT start.hdf init.hdf refine --ou=120 --rs=1 --xr='16 8 4' --ts='4 2 1' --delta='4 3 2' --an='-1 30 10' --snr=0.08 --maxit=2 --ref_a=S --cutoff=1 --MPI --olmask=80 --lmask=280 --protos='1' --hpars="-166.5 27.6" --hsearch='0.0 50.0' --oplane=10 --recon_pad=2 --full_output --sort > refine.log
