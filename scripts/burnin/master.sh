#!/bin/bash

sbatch -p csde -A csde --array=1-7 --nodes=1 --ntasks-per-node=16 --time=1:00:00 --mem=55G --job-name=s100 --export=ALL,SIMNO=100,NJOBS=7,NSIMS=100 runsim.sh
