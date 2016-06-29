#!/bin/bash

## Burnin model
qsub -q batch -t 1-10 -m n -v SIMNO=1002 runsim.burn.sh
