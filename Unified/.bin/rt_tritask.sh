#!/bin/sh

#module avail

module load compilers/cplusplus/gnu/4.4.6
module load parallel/mpi.mvapich2/1.7-r5225
module load launcher/slurm

module list

#SBATCH --time=1439
#SBATCH --partition=work
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=16

srun ./tritask > log.txt