#!/bin/bash
export OMP_PROC_BIND=TRUE
g++ -std=c++11 -Ofast  -fopenmp -mfma  -ftree-vectorize -ffast-math  -funroll-all-loops -mavx2 -mtune=native -march=native 2-D.c -o answer
