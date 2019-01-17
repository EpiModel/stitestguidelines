#!/bin/bash

cd scripts/estimation/abc
scp -r * hyak:/gscratch/csde/sjenness/stitnt
scp hyak:/gscratch/csde/sjenness/epiabc/data/*.* data/


sbatch -p csde -A csde --array=1 --export=ALL,wave=0 runsim.sh

sbatch -p csde -A csde --array=1-32 --export=ALL,wave=0 runsim.sh


sbatch -p csde -A csde --array=1-25 --export=ALL,wave=1 runsim.sh
sbatch -p csde -A csde --array=1-25 --export=ALL,wave=2 runsim.sh
sbatch -p csde -A csde --array=1-25 --export=ALL,wave=3 runsim.sh
sbatch -p csde -A csde --array=1-25 --export=ALL,wave=4 runsim.sh
sbatch -p csde -A csde --array=1-25 --export=ALL,wave=5 runsim.sh
