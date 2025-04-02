#!/bin/bash
#SBATCH -JSi8sternTol
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 10 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t24:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

# srun ../../../lib/sparc -name Si8 > out.log
srun ../../lib/rpacalc -name Si8 > out.log
mv Si8.out Si8_sternTol1e-2.out

sed -i 's/TOL_STERN_RES: 1e-2/TOL_STERN_RES: 3e-2/' Si8.rpa

srun ../../lib/rpacalc -name Si8 > out.log
mv Si8.out Si8_sternTol3e-2.out

sed -i 's/TOL_STERN_RES: 3e-2/TOL_STERN_RES: 1e-1/' Si8.rpa

srun ../../lib/rpacalc -name Si8 > out.log
mv Si8.out Si8_sternTol1e-1.out

sed -i 's/TOL_STERN_RES: 1e-1/TOL_STERN_RES: 2e-1/' Si8.rpa

srun ../../lib/rpacalc -name Si8 > out.log
mv Si8.out Si8_sternTol2e-1.out

