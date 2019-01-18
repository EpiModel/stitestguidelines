#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=1:00:00
#SBATCH --mem=100G

source ~/loadR.sh
Rscript sim.R
