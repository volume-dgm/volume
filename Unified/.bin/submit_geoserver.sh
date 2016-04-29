#!/bin/bash
#PBS -l walltime=36:00:00,nodes=1:ppn=12
#PBS -N my_job
#PBS -q batch

cd $PBS_O_WORKDIR
OMP_NUM_THREADS=$PBS_NUM_PPN ./task > task.log