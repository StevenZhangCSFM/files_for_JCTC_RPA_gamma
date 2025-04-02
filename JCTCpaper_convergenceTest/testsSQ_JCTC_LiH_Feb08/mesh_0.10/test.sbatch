#!/bin/bash
#SBATCH -JLiHmesh0.10
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 10 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t10:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

srun ../../lib/rpacalc -name LiH > out.log