#!/bin/bash
#SBATCH -JC2165Ha
#SBATCH --account=hive-phanish6
#SBATCH -N 3 --ntasks-per-node=3
#SBATCH --mem-per-cpu=10G
#SBATCH -t1:00:00
#SBATCH -phive-interact
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

srun ../../abinit-9.8.3_forModify/tmp/src/98_main/abinit< C2.files