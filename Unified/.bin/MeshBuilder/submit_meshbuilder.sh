#!/bin/sh

# It`s script for novisibirsk cluster execution

#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=16000m

cd $PBS_O_WORKDIR
./meshbuilder

cd meshes/
cp *.params ../../Task/meshes/
cp *.mesh ../../Task/meshes/
cd ..