#!/bin/bash
#SBATCH -JAbiC2120Ha
#SBATCH --account=hive-phanish6
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH --mem-per-cpu=12G
#SBATCH -t3:00:00
#SBATCH -phive
#SBATCH -oReport-%j.out

cd $SLURM_SUBMIT_DIR

module purge
module load xalt
module load intel
module load mvapich2
module load fftw
module load hdf5
module load netcdf-fortran
module load netcdf-c

srun ../../abinit-9.8.3_original/tmp/src/98_main/abinit < C2.files