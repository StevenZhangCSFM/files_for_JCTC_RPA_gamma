#!/bin/bash
#SBATCH -JC8mesh0.32
#SBATCH --account=hive-phanish6
#SBATCH -N 5 --ntasks-per-node=24
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

srun ../../lib/rpacalc -name C8 > out.log