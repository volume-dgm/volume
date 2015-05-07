#!/bin/sh

# It`s script for novisibirsk cluster execution

#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:ompthreads=1

cd $PBS_O_WORKDIR
./resultcombiner