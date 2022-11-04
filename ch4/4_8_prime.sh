#!/bin/sh
#PBS -P ACD110146
#PBS -N prime
#PBS -q ctest
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l place=scatter
#PBS -l walltime=00:01:00
#PBS -j n
module purge
module load intel/2018_u1
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
date

for ((np=1; np<=8; np=np+1)); do
  echo -n "[$np]"
  mpirun -np $np ./prime
done

mpirun -np 1 ./prime
