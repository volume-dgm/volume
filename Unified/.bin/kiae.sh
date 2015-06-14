#!/bin/sh
#SBATCH -D .
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -t 01:00:00
#SBATCH -p hpc2-16g-3d

`pwd`/meshbuilder | tee meshbuilder.log."$SLURM_JOBID"

cd meshes/
cp *.params ../../Task/meshes/
cp *.mesh ../../Task/meshes/
cd ..
