#!/bin/bash

qsub -q batch -t 1-4 -m n -v SIMNO=1000,NJOBS=4,COV=0,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1001,NJOBS=4,COV=0.5,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1002,NJOBS=4,COV=0.9,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1003,NJOBS=4,COV=0,PSTIINT=182,RC=0.5 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1004,NJOBS=4,COV=0.5,PSTIINT=182,RC=0.5 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1005,NJOBS=4,COV=0.9,PSTIINT=182,RC=0.5 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1006,NJOBS=4,COV=0,PSTIINT=182,RC=1 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1007,NJOBS=4,COV=0.5,PSTIINT=182,RC=1 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1008,NJOBS=4,COV=0.9,PSTIINT=182,RC=1 runsim.fu.sh
