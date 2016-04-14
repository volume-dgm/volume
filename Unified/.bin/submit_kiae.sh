#!/bin/sh
#SBATCH -D .
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -t 72:00:00
#SBATCH -p hpc2-16g-3d 
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive

export OMP_NUM_THREADS=8

module load openmpi/latest

$MPIRUN `pwd`/task 
