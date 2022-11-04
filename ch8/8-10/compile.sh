#!/bin/bash
mpicc -o code ./code.c -lm -std=c11
gcc -o validate ./validate.c -std=c11 -lm
gcc -o gen-vector ./gen-vector.c -std=c11 -lm
gcc -o gen-double-matrix ./gen-double-matrix.c -std=c11 -lm
