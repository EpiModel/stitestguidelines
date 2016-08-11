#!/bin/bash

qsub -q batch -t 1-2 -m n -v SIMNO=1000,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1001,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.1 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1002,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.2 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1003,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.3 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1004,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.4 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1005,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.5 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1006,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.6 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1007,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.7 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1008,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.8 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1009,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=0.9 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1010,NJOBS=2,COV=0,PSTIINT=182,RC=0,ASYMPT=1 runsim.fu.sh
