#!/bin/bash
#SBATCH -JSi2KS
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 1 --ntasks-per-node=12
#SBATCH --mem-per-cpu=7G
#SBATCH -t1:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

srun ../../../lib/sparc -name Si2 > out.log