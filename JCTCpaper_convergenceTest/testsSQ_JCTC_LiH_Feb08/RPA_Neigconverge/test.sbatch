#!/bin/bash
#SBATCH -JLiHNeig
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 10 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t6:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

# srun ../../lib/rpacalc -name LiH > out.log
# mv LiH.out LiH_Neig320.out

# sed -i 's/N_NUCHI_EIGS: 320 # for testing /N_NUCHI_EIGS: 640 # for testing /' LiH.rpa

# srun ../../lib/rpacalc -name LiH > out.log
# mv LiH.out LiH_Neig640.out

# sed -i 's/N_NUCHI_EIGS: 640 # for testing /N_NUCHI_EIGS: 1280 # for testing /' LiH.rpa

# srun ../../lib/rpacalc -name LiH > out.log
# mv LiH.out LiH_Neig1280.out

# sed -i 's/N_NUCHI_EIGS: 1280 # for testing /N_NUCHI_EIGS: 2240 # for testing /' LiH.rpa

# srun ../../lib/rpacalc -name LiH > out.log
# mv LiH.out LiH_Neig2240.out

# sed -i 's/N_NUCHI_EIGS: 2240 # for testing /N_NUCHI_EIGS: 3200 # for testing /' LiH.rpa

# srun ../../lib/rpacalc -name LiH > out.log
# mv LiH.out LiH_Neig3200.out

# sed -i 's/N_NUCHI_EIGS: 3200 # for testing /N_NUCHI_EIGS: 4160 # for testing /' LiH.rpa

srun ../../lib/rpacalc -name LiH > out.log
mv LiH.out LiH_Neig4800.out
sed -i 's/N_NUCHI_EIGS: 4800 # for testing /N_NUCHI_EIGS: 5600 # for testing /' LiH.rpa
srun ../../lib/rpacalc -name LiH > out.log
mv LiH.out LiH_Neig5600.out