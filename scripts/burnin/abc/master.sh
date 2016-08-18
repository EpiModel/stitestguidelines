#!/bin/bash

qsub -q batch -m ae -v NSIM=100,PACC=0.01 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.01 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=500,PACC=0.01 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.02 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.02 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.05 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.10 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=200,PACC=0.10 runsim.burn.abcsmc.sh
