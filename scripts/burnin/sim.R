
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

                   rgc.tprob = p[1]+0.00445,
                   ugc.tprob = p[2]+0.00445,
                   rct.tprob = p[3]-0.0006,
                   uct.tprob = p[4]-0.0006,

                   rgc.asympt.rate = 1/(24.78753*7),
                   ugc.asympt.rate = 1/(24.78753*7),
                   rct.asympt.rate = 1/(44.28232*7),
                   uct.asympt.rate = 1/(44.28232*7),

                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.80,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.48,

                   hiv.rgc.rr = p[5]+0.02,
                   hiv.ugc.rr = p[6]+0.02,
                   hiv.rct.rr = p[5]+0.02,
                   hiv.uct.rr = p[6]+0.02,

                   ai.scale = p[7]+0.0032,
                   ai.scale.pospos = p[7]+0.0032,

                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.ct.hivpos.coverage = 0.61,

                   stitest.start = 1,
                   sti.correlation.time = 0)
init <- init_msm(nwstats = st,
                 prev.B = p[8],
                 prev.W = p[8],
                 prev.ugc = p[9]-0.0016,
                 prev.rgc = p[9]-0.0016,
                 prev.uct = p[10]+0.0018,
                 prev.rct = p[10]+0.0018,
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
