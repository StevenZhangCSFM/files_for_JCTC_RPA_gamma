#!/bin/bash
#SBATCH -JSi8mesh0.45
#SBATCH --account=hive-phanish6
#SBATCH -N 3 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t1:00:00
#SBATCH -phive-interact
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

# srun ../../../lib/sparc -name Si8 > out.log
srun ../../lib/rpacalc -name Si8 > out.log