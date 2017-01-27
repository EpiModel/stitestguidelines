
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
                   ai.scale = 1.11915125,
                   prep.start = 5000,
                   partlist.start = 5000,
                   stitest.start = 5000,
                   
                   hiv.rgc.rr = 2.780673,
                   hiv.ugc.rr = 1.732363,
                   hiv.rct.rr = 2.780673,
                   hiv.uct.rr = 1.732363,
                   hiv.syph.rr = 2.00892218,
                   syph.hiv.rr = 2.17038393,
                   
                   syph.tprob = 0.01983336,
                   rgc.tprob = 0.3928965, # was 0.357698
                   ugc.tprob = 0.24297633, # was 0.248095
                   rct.tprob = 0.29367628, # was 0.321597
                   uct.tprob = 0.25309465, # was 0.212965
                   
                   prep.coverage = 0,
                   ept.coverage = 0,
                   
                   rgc.sympt.prob = 0.16, # Beck
                   ugc.sympt.prob = 0.90, # Beck
                   rct.sympt.prob = 0.14, # Beck
                   uct.sympt.prob = 0.58, # Beck
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182)
                   
init <- init_msm(nwstats = st,
                 prev.B = 0.10, 
                 prev.W = 0.10,
                 prev.ugc = 0.013,
                 prev.rgc = 0.013,
                 prev.uct = 0.013,
                 prev.rct = 0.013,
                 prev.syph.B = 0.01,
                 prev.syph.W = 0.01)
control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
