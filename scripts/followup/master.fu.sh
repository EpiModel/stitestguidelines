#!/bin/bash

qsub -q batch -t 1-2 -m n -v SIMNO=1000,NJOBS=2,COV=0,PSTIINT=7,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1001,NJOBS=2,COV=0,PSTIINT=21,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1002,NJOBS=2,COV=0,PSTIINT=35,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1003,NJOBS=2,COV=0,PSTIINT=49,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1004,NJOBS=2,COV=0,PSTIINT=63,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1005,NJOBS=2,COV=0,PSTIINT=77,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1006,NJOBS=2,COV=0,PSTIINT=91,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1007,NJOBS=2,COV=0,PSTIINT=105,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1008,NJOBS=2,COV=0,PSTIINT=119,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1009,NJOBS=2,COV=0,PSTIINT=133,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1010,NJOBS=2,COV=0,PSTIINT=147,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1011,NJOBS=2,COV=0,PSTIINT=161,RC=0,ASYMPT=0 runsim.fu.sh
qsub -q batch -t 1-2 -m n -v SIMNO=1012,NJOBS=2,COV=0,PSTIINT=175,RC=0,ASYMPT=0 runsim.fu.sh
