#!/bin/bash
#SBATCH -JLiH442RPA
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 24 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t30:00:00
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