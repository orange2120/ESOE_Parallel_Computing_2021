#!/bin/bash
#PBS -P ACD110146
#PBS -N code
#PBS -l select=1:ncpus=20:mpiprocs=20
#PBS -l walltime=00:30:00
#PBS -q ctest
#PBS -j n
module purge
module load intel/2018_u1
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

date
for ((np=1; np<=20; np=np+1)); do
  mpirun -n $np ./code 10000000
done
