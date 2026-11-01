#!/bin/bash
#SBATCH -JeigLiH111
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem-per-cpu=7G
#SBATCH -t01:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

module reset

module load intel-oneapi-mkl
module load intel-oneapi-mpi
# module load mvapich2

srun ../../../../lib/sparc -name LiH > out.log
