#!/bin/bash
#SBATCH --job-name="myo6actin"


cd $SLURM_SUBMIT_DIR

source /data/Alushinlab/environment/.bashrc
#source /data/Alushinlab/environment/EMAN1-1.9/eman.bashrc

#load required module, cross fingers
module load EMAN2/2.1

#execute

mpirun -np 64 /data/Alushinlab/lab_scripts/EMAN2/EMAN2_actin2/refine.py start.hdf init.hdf refine --ou=120 --rs=1 --xr='16 8 4' --ts='4 2 1' --delta='4 3 2' --an='-1 30 10' --snr=0.08 --maxit=2 --ref_a=S --cutoff=1.5 --MPI --olmask=80 --lmask=280 --protos='1' --hpars='-166.5 27.6' --hsearch='0.0 50.0' --oplane=10 --recon_pad=2 --full_output --sort > refine.log

