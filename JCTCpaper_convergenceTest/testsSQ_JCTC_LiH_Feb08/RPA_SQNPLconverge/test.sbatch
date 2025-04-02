#!/bin/bash
#SBATCH -JLiHSQ
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 10 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t12:00:00
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
mv LiH.out LiH_SQ5.out

sed -i 's/SQ_NPL_RPA: 5/SQ_NPL_RPA: 4/' LiH.rpa

srun ../../lib/rpacalc -name LiH > out.log
mv LiH.out LiH_SQ4.out

sed -i 's/SQ_NPL_RPA: 4/SQ_NPL_RPA: 3/' LiH.rpa

srun ../../lib/rpacalc -name LiH > out.log
mv LiH.out LiH_SQ3.out

sed -i 's/SQ_NPL_RPA: 3/SQ_NPL_RPA: 2/' LiH.rpa

srun ../../lib/rpacalc -name LiH > out.log
mv LiH.out LiH_SQ2.out