#!/bin/bash

sbatch -p csde -A csde --array=1-32 --nodes=1 --ntasks-per-node=16 --time=2:00:00 --mem=55G --job-name=s3000 --export=ALL,SIMNO=3000,NJOBS=32,NSIMS=500 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-32 --nodes=1 --ntasks-per-node=16 --time=2:00:00 --mem=55G --job-name=s3001 --export=ALL,SIMNO=3001,NJOBS=32,NSIMS=500 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-63 --nodes=1 --ntasks-per-node=16 --time=2:00:00 --mem=55G --job-name=s3002 --export=ALL,SIMNO=3002,NJOBS=63,NSIMS=1000 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-125 --nodes=1 --ntasks-per-node=16 --time=2:00:00 --mem=55G --job-name=s3003 --export=ALL,SIMNO=3003,NJOBS=125,NSIMS=2000 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-125 --nodes=1 --ntasks-per-node=16 --time=2:00:00 --mem=55G --job-name=s3004 --export=ALL,SIMNO=3004,NJOBS=125,NSIMS=2000 runsim.sh

sbatch -p ckpt -A csde-ckpt --array=1-250 --nodes=1 --ntasks-per-node=16 --time=2:00:00 --mem=55G --job-name=s3005 --export=ALL,SIMNO=3005,NJOBS=250,NSIMS=4000 runsim.sh
