
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiABC"))

## Environmental Arguments
pull_env_vars()

## Parameters
load("est/nwstats.20k.rda")
p <- get_posterior(wave = 27, input = "est/")
p <- colMeans(p$param)


param <- param_msm(nwstats = st,

                   rgc.tprob = 0.5364416,
                   ugc.tprob = 0.434692,
                   rct.tprob = 0.2493814,
                   uct.tprob = 0.1944415,

                   rgc.asympt.rate = 1/(24.78753*7),
                   ugc.asympt.rate = 1/(24.78753*7),
                   rct.asympt.rate = 1/(44.28232*7),
                   uct.asympt.rate = 1/(44.28232*7),

                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.80,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.48,

                   hiv.rgc.rr = 2.175918,
                   hiv.ugc.rr = 1.564797,
                   hiv.rct.rr = 2.175918,
                   hiv.uct.rr = 1.564797,

                   ai.scale = 1.061338,
                   ai.scale.pospos = 1.061338,

                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.ct.hivpos.coverage = 0.61,

                   stitest.start = 1,
                   sti.correlation.time = 0)
init <- init_msm(nwstats = st,
                 prev.B = 0.149006,
                 prev.W = 0.149006,
                 prev.ugc = 0.001471584,
                 prev.rgc = 0.001471584,
                 prev.uct = 0.007572175,
                 prev.rct = 0.007572175,
                 prev.syph.B = 0,
                 prev.syph.W = 0)
control <- control_msm(simno = fsimno,
                       nsteps = 52*50,
                       nsims = ncores,
                       ncores = ncores,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.20k.rda", param, init, control,
           save.min = FALSE, save.max = TRUE,
           compress = TRUE, verbose = FALSE)

# # Post-Merging
# vars <- c("ir100.gc", "ir100.ct", "i.prev")
# process_simfiles(simno = simno, min.n = njobs, nsims = nsims,
#                  vars = vars, delete.sub = TRUE)
