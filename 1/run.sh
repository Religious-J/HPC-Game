#!/bin/bash
#SBATCH --partition=compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1

srun ./program > output.dat
