#!/bin/bash
#SBATCH -JeigLiH111
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 3 --ntasks-per-node=12
#SBATCH --mem-per-cpu=7G
#SBATCH -t00:10:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

module reset

module load intel-oneapi-mkl
module load intel-oneapi-mpi

srun ../../lib/sparc -name LiH > out.log

# # doing data analysis, finding the time spent on each part
# module load anaconda3 # for using numpy

# python ../../scriptSQ.py > summary.log

# python ../../statistic.py 36 1 1 > summary2.log