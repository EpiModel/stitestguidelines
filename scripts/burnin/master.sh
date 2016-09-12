#!/bin/bash

# runs abc fit
qsub -q batch runsim.burn.abcsmc2.sh


# runs burnin model
qsub -q bf -t 1-16 -v SIMNO=100,NJOBS=16 runsim.burn.sh
qsub -q bf -t 17-32 -v SIMNO=100,NJOBS=16 runsim.burn.sh
