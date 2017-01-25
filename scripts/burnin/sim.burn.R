
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
                   hiv.uct.rr = 1.732363,
                   hiv.syph.rr = 2.00,
                   syph.hiv.rr = 2.40,
                   

                   syph.tprob = 0.020,
                   rgc.tprob = 0.40, # was 0.357698
                   ugc.tprob = 0.30, # was 0.248095
                   rct.tprob = 0.25, # was 0.321597
                   uct.tprob = 0.15, # was 0.212965

                   ai.scale = 1.12,
                   prep.coverage = 0,
                   ept.coverage = 0,
                   
                   rgc.sympt.prob = 0.16, # Beck
                   ugc.sympt.prob = 0.90, # Beck
                   rct.sympt.prob = 0.14, # Beck
                   uct.sympt.prob = 0.58, # Beck
                   
                   stitest.int = 182) # adjustable for 3 or 6 months)
init <- init_msm(nwstats = st)
control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
