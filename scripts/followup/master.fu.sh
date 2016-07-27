#!/bin/bash

qsub -q batch -t 1-4 -m n -v SIMNO=1000,COV=0,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1001,COV=0.5,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1002,COV=0.9,PSTIINT=182,RC=0 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1003,COV=0,PSTIINT=182,RC=0.5 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1004,COV=0.5,PSTIINT=182,RC=0.5 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1005,COV=0.9,PSTIINT=182,RC=0.5 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1006,COV=0,PSTIINT=182,RC=1 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1007,COV=0.5,PSTIINT=182,RC=1 runsim.fu.sh
qsub -q batch -t 1-4 -m n -v SIMNO=1008,COV=0.9,PSTIINT=182,RC=1 runsim.fu.sh
