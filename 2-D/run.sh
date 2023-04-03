#!/bin/bash
srun -N 1 -n 1 -c 1 -p compute  ./program $1 >> output.dat
seff "$(cat job_id.dat)" > seff.dat