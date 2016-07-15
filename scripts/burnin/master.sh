#!/bin/bash

qsub -q batch -t 1-12 runsim.burn.abcr.sh

qsub -q bf -t 1-200 runsim.burn.abcr.sh
qsub -q bf -t 201-500 runsim.burn.abcr.sh
qsub -q bf -t 501-600 runsim.burn.abcr.sh
