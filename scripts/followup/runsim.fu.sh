#!/bin/bash

### User specs
#PBS -N followup
#PBS -l nodes=1:ppn=16,mem=50gb,feature=16core,walltime=01:00:00
#PBS -o /gscratch/csde/kweiss2/sti/out
#PBS -e /gscratch/csde/kweiss2/sti/out
#PBS -j oe
#PBS -d /gscratch/csde/kweiss2/sti
#PBS -m n

### Standard specs
HYAK_NPE=$(wc -l < $PBS_NODEFILE)
HYAK_NNODES=$(uniq $PBS_NODEFILE | wc -l )
HYAK_TPN=$((HYAK_NPE/HYAK_NNODES))
NODEMEM=`grep MemTotal /proc/meminfo | awk '{print $2}'`
NODEFREE=$((NODEMEM-2097152))
MEMPERTASK=$((NODEFREE/HYAK_TPN))
ulimit -v $MEMPERTASK
export MX_RCACHE=0

### Modules
module load msropen_3.3.2

### App
Rscript sim.fu.R
