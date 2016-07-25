#!/bin/bash

qsub -q batch -t 1-25 -v BATCH=2,SETSIZE=100 runsim.burn.abcsmc.sh
qsub -q bf -t 26-100 -v BATCH=2,SETSIZE=100 runsim.burn.abcsmc.sh

qsub -q batch -t 1-25 -v BATCH=3,SETSIZE=100 runsim.burn.abcsmc.sh
qsub -q bf -t 26-100 -v BATCH=3,SETSIZE=100 runsim.burn.abcsmc.sh

qsub -q batch -t 1-25 -v BATCH=4,SETSIZE=100 runsim.burn.abcsmc.sh
qsub -q bf -t 26-100 -v BATCH=4,SETSIZE=100 runsim.burn.abcsmc.sh
qsub -q bf -t 28 -v BATCH=4,SETSIZE=100 runsim.burn.abcsmc.sh

qsub -q batch -t 1-25 -v BATCH=5,SETSIZE=100 runsim.burn.abcsmc.sh
qsub -q bf -t 26-100 -v BATCH=5,SETSIZE=100 runsim.burn.abcsmc.sh

qsub -q batch runsim.burn.abcsmc2.sh

qsub -q batch -t 1-50 -v BATCH=6,SETSIZE=100 runsim.burn.abcsmc.sh
qsub -q bf -t 51-250 -v BATCH=6,SETSIZE=100 runsim.burn.abcsmc.sh

qsub -q batch -t 1-7 -v SIMNO=100 runsim.burn.sh
