#!/bin/bash

qsub -q batch -t 1-7 -m n -v SIMNO=1000,NJOBS=7,COV=0.4,PSTIINT=7,RC=0 runsim.fu.sh
qsub -q batch -t 1-7 -m n -v SIMNO=1001,NJOBS=7,COV=0.4,PSTIINT=14,RC=0 runsim.fu.sh
qsub -q batch -t 1-7 -m n -v SIMNO=1002,NJOBS=7,COV=0.4,PSTIINT=28,RC=0 runsim.fu.sh
qsub -q batch -t 1-7 -m n -v SIMNO=1003,NJOBS=7,COV=0.4,PSTIINT=56,RC=0 runsim.fu.sh
qsub -q batch -t 1-7 -m n -v SIMNO=1004,NJOBS=7,COV=0.4,PSTIINT=91,RC=0 runsim.fu.sh
qsub -q batch -t 1-7 -m n -v SIMNO=1005,NJOBS=7,COV=0.4,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-7 -m n -v SIMNO=1006,NJOBS=7,COV=0.4,PSTIINT=364,RC=0 runsim.fu.sh
