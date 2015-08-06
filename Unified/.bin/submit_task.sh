#!/bin/sh

#It`s for novosibirsk cluster execution

#PBS -l walltime=48:00:00
#PBS -l select=16:ncpus=12:mpiprocs=1:ompthreads=12:mem=2000m

cd $PBS_O_WORKDIR
mpirun -hostfile $PBS_NODEFILE --bind-to-core ./task > log.txt

