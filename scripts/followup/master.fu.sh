#!/bin/bash

qsub -q batch -t 1-4 -m n -v SIMNO=1000,NJOBS=4,COV=0,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1001,NJOBS=4,COV=0.4,PSTIINT=182,RC=0 runsim.fu.sh

