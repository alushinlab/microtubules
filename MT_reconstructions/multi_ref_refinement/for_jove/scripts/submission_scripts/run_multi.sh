#!/bin/bash
#PBS -N myo7


cd $PBS_O_WORKDIR

source /data/Alushinlab/environment/.bashrc
source /usr/local/EMAN1-1.9/eman.bashrc

#load required module, cross fingers
module load eman2

#execute

/usr/local/apps/EMAN2/openmpi/1.6.5/gnu/eth/bin/mpirun -machinefile $PBS_NODEFILE -np 32 /data/Alushinlab/lab_scripts/EMAN2/EMAN2_actin2/refine.py start.hdf init.hdf refine --ou=120 --rs=1 --xr='16 8 4' --ts='4 2 1' --delta='4 3 2' --an='-1 30 10' --snr=0.08 --maxit=2 --ref_a=S --cutoff=1 --MPI --olmask=80 --lmask=280 --protos='1 1' --hpars='-166.5 27.6 -166.5 27.6' --hsearch='0.0 50.0' --oplane=10 --recon_pad=2 --full_output --sort > refine.log
