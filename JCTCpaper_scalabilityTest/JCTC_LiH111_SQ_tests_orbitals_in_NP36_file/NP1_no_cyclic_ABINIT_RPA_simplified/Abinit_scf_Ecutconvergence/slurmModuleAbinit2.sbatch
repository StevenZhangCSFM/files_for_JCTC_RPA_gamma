#!/bin/bash
#SBATCH -JAbinitLiH
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=40G
#SBATCH -t01:00:00
#SBATCH -oReport-%j.out

cd $SLURM_SUBMIT_DIR

module load gcc/12.3.0
module load openmpi/4.1.5
module load hdf5/1.8.22-ompi
module load netcdf-c/4.9.2-ompi-hdf5-1.8
module load intel-oneapi-mkl
module load fftw/3.3.10
module load c-blosc/1.21.2-abinit
module load abinit

srun abinit < LiH.files