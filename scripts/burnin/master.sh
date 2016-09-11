#!/bin/bash

# runs abc fit
qsub -q batch runsim.burn.abcsmc2.sh


# runs burnin model
qsub -q batch -t 1-32 -v SIMNO=100,NJOBS=32 runsim.burn.sh
