#!/bin/sh
#PBS -P ACD110146
#PBS -N pi_calc
#PBS -q ctest
#PBS -l select=1:ncpus=24:mpiprocs=24
#PBS -l place=scatter
#PBS -l walltime=00:01:00
#PBS -j n
module purge
module load intel/2018_u1
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR
date

for ((np=1; np<=24; np=np+1)); do
  echo "Proc [$np]"
  mpirun -np $np ./pi_calc
done

mpirun ./pi_calc
