
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
                   prep.start = 5000,
                   rgc.dur.asympt = 35.11851,
                   ugc.dur.asympt = 35.11851,
                   rct.dur.asympt = 44.24538,
                   uct.dur.asympt = 44.24538,

                   hiv.rgc.rr = 2.780673,
                   hiv.ugc.rr = 1.732363,
                   hiv.rct.rr = 2.780673,
                   hiv.uct.rr = 1.732363)
init <- init_msm(nwstats = st)
control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
