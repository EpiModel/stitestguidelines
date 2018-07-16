#!/bin/bash

## New Slurm -------------------------------------------------------------------
sbatch -p csde -A csde --export=NSIM=100,PACC=0.01 runsim.sh
sbatch -p csde -A csde --export=NSIM=100,PACC=0.05 runsim.sh
sbatch -p csde -A csde --export=NSIM=100,PACC=0.10 runsim.sh

sbatch -p ckpt  -A csde-ckpt --export=NSIM=100,PACC=0.01 runsim.sh
sbatch -p ckpt  -A csdee-ckpt --export=NSIM=100,PACC=0.05 runsim.sh
sbatch -p ckpt  -A csdee-ckpt --export=NSIM=100,PACC=0.10 runsim.sh


## Old Hyak --------------------------------------------------------------------

qsub -q batch -m ae -v NSIM=100,PACC=0.01 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.10 runsim.burn.abcsmc.sh


qsub -q batch -m ae -v NSIM=100,PACC=0.01 runsim.burn.abcsmc2.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc2.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.10 runsim.burn.abcsmc2.sh


qsub -q batch -m ae -v NSIM=100,PACC=0.01 runsim.burn.abcsmc3.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc3.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.10 runsim.burn.abcsmc3.sh

qsub -q batch -m ae -v NSIM=100,PACC=0.01 runsim.burn.abcsmc4.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc4.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.10 runsim.burn.abcsmc4.sh

qsub -q batch -m ae -v NSIM=100,PACC=0.01 runsim.burn.abcsmc4.syph.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc4.syph.sh
qsub -q batch -m ae -v NSIM=100,PACC=0.10 runsim.burn.abcsmc4.syph.sh

qsub -q int -m ae -v NSIM=100,PACC=0.05 runsim.burn.abcsmc4.syph.sh
qsub -q int -m ae -v NSIM=200,PACC=0.01 runsim.burn.abcsmc4.syph.sh
qsub -q int -m ae -v NSIM=500,PACC=0.01 runsim.burn.abcsmc4.syph.sh
