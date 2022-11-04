#!/bin/bash
#PBS -P ACD110146
#PBS -N area
#PBS -q ctest
#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l place=scatter
#PBS -l walltime=00:03:00
#PBS -j n

n_sample=1000000

module purge
module load intel/2018_u1
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
date

for ((np=1; np<=24; np=np+1)); do
	#   echo "Proc [$np]"
	  mpirun -np $np ./code $n_sample;
done
