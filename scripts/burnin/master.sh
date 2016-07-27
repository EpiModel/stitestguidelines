#!/bin/bash

# runs abc fit
qsub -q batch runsim.burn.abcsmc2.sh


# runs burnin model
qsub -q batch -t 1-7 -v SIMNO=300 runsim.burn.sh



qsub -q batch runsim.burn.abcsmc3.sh
