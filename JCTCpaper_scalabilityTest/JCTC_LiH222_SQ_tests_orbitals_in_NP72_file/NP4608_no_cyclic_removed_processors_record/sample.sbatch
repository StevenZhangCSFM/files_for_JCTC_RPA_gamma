#!/bin/bash
#SBATCH -JeigLiH222
#SBATCH --account=hive-phanish6
#SBATCH -N 192 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t01:00:00
#SBATCH -phive
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

module reset

module load intel-oneapi-mkl
module load intel-oneapi-mpi

srun ../../lib/rpacalc -name LiH > out.log

# doing data analysis, finding the time spent on each part
module load anaconda3 # for using numpy

python ../../scriptSQ.py > summary.log

python ../../statistic.py 4608 1 1 > summary2.log