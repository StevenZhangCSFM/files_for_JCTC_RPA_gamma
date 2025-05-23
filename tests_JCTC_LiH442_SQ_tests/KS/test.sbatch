#!/bin/bash
#SBATCH -JLiH442KS
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 2 --ntasks-per-node=12
#SBATCH --mem-per-cpu=7G
#SBATCH -t00:30:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

srun ../../lib/sparc -name LiH > out.log