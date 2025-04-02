#!/bin/bash
#SBATCH -JC8Neig
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 10 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t8:00:00
#SBATCH -oReport-%j.out
cd $SLURM_SUBMIT_DIR

# module load intel
# module load mvapich2

# module load gcc
# module load mkl
# module load mvapich2

module load intel-oneapi-mkl
module load intel-oneapi-mpi

# srun ../../lib/rpacalc -name C8 > out.log
# mv C8.out C8_Neig640.out

# sed -i 's/N_NUCHI_EIGS: 640 # for testing /N_NUCHI_EIGS: 1280 # for testing /' C8.rpa

# srun ../../lib/rpacalc -name C8 > out.log
# mv C8.out C8_Neig1280.out

# sed -i 's/N_NUCHI_EIGS: 1280 # for testing /N_NUCHI_EIGS: 2560 # for testing /' C8.rpa

# srun ../../lib/rpacalc -name C8 > out.log
# mv C8.out C8_Neig2560.out

# sed -i 's/N_NUCHI_EIGS: 2560 # for testing /N_NUCHI_EIGS: 4480 # for testing /' C8.rpa

# srun ../../lib/rpacalc -name C8 > out.log
# mv C8.out C8_Neig4480.out

# sed -i 's/N_NUCHI_EIGS: 4480 # for testing /N_NUCHI_EIGS: 6400 # for testing /' C8.rpa

# srun ../../lib/rpacalc -name C8 > out.log
# mv C8.out C8_Neig6400.out

# # sed -i 's/N_NUCHI_EIGS: 6400 # for testing /N_NUCHI_EIGS: 8320 # for testing /' C8.rpa

srun ../../lib/rpacalc -name C8 > out.log
mv C8.out C8_Neig9600.out
sed -i 's/N_NUCHI_EIGS: 9600 # for testing /N_NUCHI_EIGS: 11200 # for testing /' C8.rpa
srun ../../lib/rpacalc -name C8 > out.log
mv C8.out C8_Neig11200.out