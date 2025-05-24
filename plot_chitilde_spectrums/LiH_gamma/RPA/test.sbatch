#!/bin/bash
#SBATCH -JLiH0.15
#SBATCH --account=hive-phanish6
#SBATCH -N 2 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t00:30:00
#SBATCH -phive-interact
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

# srun ../../../lib/sparc -name LiH > out.log
srun ../../lib/rpacalc -name LiH > out.log