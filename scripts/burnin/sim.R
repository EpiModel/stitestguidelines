
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
pull_env_vars()

## Parameters
load("est/nwstats.rda")
param <- param_msm(nwstats = st,

                   rgc.tprob = 0.570, # x2
                   ugc.tprob = 0.500, # x3
                   rct.tprob = 0.245, # x4
                   uct.tprob = 0.205, # x5

                   rgc.asympt.int = 21*7, # x6
                   ugc.asympt.int = 21*7, # x6
                   rct.asympt.int = 46,
                   uct.asympt.int = 46,

                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.80,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.48,

                   hiv.rgc.rr = 1.97,
                   hiv.ugc.rr = 1.48,
                   hiv.rct.rr = 1.97,
                   hiv.uct.rr = 1.48,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.ct.hivpos.coverage = 0.61,

                   stitest.start = 1,
                   sti.correlation.time = 0)
init <- init_msm(nwstats = st,
                 prev.B = 0.14,
                 prev.W = 0.14,
                 prev.ugc = 0.005,
                 prev.rgc = 0.005,
                 prev.uct = 0.015,
                 prev.rct = 0.015,
                 prev.syph.B = 0,
                 prev.syph.W = 0)
control <- control_msm(simno = fsimno,
                       nsteps = 500,
                       nsims = ncores,
                       ncores = ncores,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
           save.min = TRUE, save.max = FALSE,
           compress = FALSE, verbose = FALSE)

# Post-Merging
vars <- c("ir100.gc", "ir100.ct", "i.prev", "num")
process_simfiles(simno = simno, min.n = njobs, nsims = nsims,
                 vars = vars, delete.sub = TRUE)
