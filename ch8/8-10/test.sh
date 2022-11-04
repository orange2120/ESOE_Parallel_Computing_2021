# Usage: ./test.sh <the size of the square matrix>
sh generate.sh $1
mpirun -np 25 ./code matrix vector output
./validate matrix vector output
