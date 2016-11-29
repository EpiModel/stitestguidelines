
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
load("est/abc.syph.parms.rda")
for (i in seq_along(mean.p)) {
    assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,
                   syph.tprob = 0.024,
                   hiv.syph.rr = 2.4,
                   syph.hiv.rr = 2.3,
                   prep.start = 5000)
init <- init_msm(nwstats = st)
control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
