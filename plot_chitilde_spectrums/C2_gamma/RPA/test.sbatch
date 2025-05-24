#!/bin/bash
#SBATCH -JC2test
#SBATCH --account=hive-phanish6
#SBATCH -N 2 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t01:00:00
#SBATCH -phive-interact
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

module reset
module load intel-oneapi-mkl
module load intel-oneapi-mpi

# srun ../../../lib/sparc -name BaTiO3_kpt
srun ../../lib/rpacalc -name C2 > out.log