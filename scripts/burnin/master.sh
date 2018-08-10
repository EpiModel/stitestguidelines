#!/bin/bash

# Newer Slurm

sbatch -p csde -A csde --array=1-16 --export=SIMNO=100,NJOBS=16 runsim.burn.sh
sbatch -p csde -A csde --array=1-16 --export=SIMNO=101,NJOBS=16 runsim.burn.sh

sbatch -p ckpt -A csde-ckpt --array=1-16 --export=SIMNO=100,NJOBS=16 runsim.burn.sh

## Old Hyak --------------------------------------------------------------------

# runs abc fit
qsub -q batch runsim.burn.abcsmc2.sh

# runs burnin model
qsub -q batch -t 1-16 -v SIMNO=100,NJOBS=16 runsim.burn.sh
qsub -q batch -t 1-16 -v SIMNO=101,NJOBS=16 runsim.burn.sh


qsub -q bf -t 1-16 -v SIMNO=100,NJOBS=16 runsim.burn.sh
qsub -q bf -t 1-32 -v SIMNO=105,NJOBS=32 runsim.burn.sh
qsub -q bf -t 1-20 -v SIMNO=201,NJOBS=20 runsim.burn.sh

