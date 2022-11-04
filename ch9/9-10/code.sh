#!/bin/bash
#PBS -P ACD110146
#PBS -N code
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l walltime=00:30:00
#PBS -q ctest
#PBS -j n
module purge
module load intel/2018_u1
cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

date
for ((np=2; np<9; np=np+1)); do
  mpirun -n $np ./code
done
mpirun -n 1 ./ans
