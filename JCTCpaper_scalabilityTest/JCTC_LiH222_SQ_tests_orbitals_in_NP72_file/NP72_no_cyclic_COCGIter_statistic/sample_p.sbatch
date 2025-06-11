#!/bin/bash
#SBATCH -JeigLiH222
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 3 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t03:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

module reset

module load intel-oneapi-mkl
module load intel-oneapi-mpi

srun ../../lib/sparc -name LiH > out.log
