#!/bin/bash
#SBATCH -JSi8Neig
#SBATCH --account=gts-phanish6-coda20
#SBATCH -N 10 --ntasks-per-node=24
#SBATCH --mem-per-cpu=7G
#SBATCH -t18:00:00
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
# srun ../../lib/rpacalc -name Si8 > out.log
# mv Si8.out Si8_Neig640.out

# sed -i 's/N_NUCHI_EIGS: 640 # for testing /N_NUCHI_EIGS: 1280 # for testing /' Si8.rpa

# srun ../../lib/rpacalc -name Si8 > out.log
# mv Si8.out Si8_Neig1280.out

# sed -i 's/N_NUCHI_EIGS: 1280 # for testing /N_NUCHI_EIGS: 2560 # for testing /' Si8.rpa

# srun ../../lib/rpacalc -name Si8 > out.log
# mv Si8.out Si8_Neig2560.out

# sed -i 's/N_NUCHI_EIGS: 2560 # for testing /N_NUCHI_EIGS: 4480 # for testing /' Si8.rpa

# srun ../../lib/rpacalc -name Si8 > out.log
# mv Si8.out Si8_Neig4480.out

# sed -i 's/N_NUCHI_EIGS: 4480 # for testing /N_NUCHI_EIGS: 6400 # for testing /' Si8.rpa

# srun ../../lib/rpacalc -name Si8 > out.log
# mv Si8.out Si8_Neig6400.out

# sed -i 's/N_NUCHI_EIGS: 6400 # for testing /N_NUCHI_EIGS: 8320 # for testing /' Si8.rpa

srun ../../lib/rpacalc -name Si8 > out.log
mv Si8.out Si8_Neig9600.out
sed -i 's/N_NUCHI_EIGS: 9600 # for testing /N_NUCHI_EIGS: 11200 # for testing /' Si8.rpa
srun ../../lib/rpacalc -name Si8 > out.log
mv Si8.out Si8_Neig11200.out