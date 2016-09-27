
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- Sys.getenv("SIMNO")
jobno <- Sys.getenv("PBS_ARRAYID")
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   prep.start = 3121)
init <- init_msm(nwstats = st)
control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
