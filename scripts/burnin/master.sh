#!/bin/bash

## Burnin model
qsub -q batch -t 1-7 -m n -v SIMNO=1000 runsim.burn.sh
