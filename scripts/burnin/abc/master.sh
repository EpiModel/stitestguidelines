#!/bin/bash

cd scripts/burnin/abc
scp -r * hyak:/gscratch/csde/sjenness/stitnt
scp hyak:/gscratch/csde/sjenness/stitnt/data/*.rda data/

sbatch -p csde -A csde --array=1-125 --export=ALL,wave=0 runsim.sh
sbatch -p csde -A csde --array=1-100 --export=ALL,wave=1 runsim.sh
sbatch -p csde -A csde --array=1-100 --export=ALL,wave=2 runsim.sh
sbatch -p csde -A csde --array=1-100 --export=ALL,wave=3 runsim.sh
sbatch -p csde -A csde --array=1-100 --export=ALL,wave=4 runsim.sh
