
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
load("est/meta.parms.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,

                   ai.scale = 1.15,

                   rgc.tprob = rgc.tprob - 0.066,
                   ugc.tprob = ugc.tprob - 0.008,
                   rct.tprob = rct.tprob - 0.09,
                   uct.tprob = uct.tprob - 0.022,

                   rgc.sympt.prob = rgc.sympt.prob,
                   ugc.sympt.prob = ugc.sympt.prob + 0.076,
                   rct.sympt.prob = rct.sympt.prob,
                   uct.sympt.prob = uct.sympt.prob,

                   rgc.dur.asympt = rgc.dur.asympt - 0.215,
                   ugc.dur.asympt = ugc.dur.asympt - 0.215,

                   rct.dur.asympt = rct.dur.asympt - 0.782,
                   uct.dur.asympt = uct.dur.asympt - 0.782,

                   gc.prob.cease = 0,
                   ct.prob.cease = 0,

                   hiv.rgc.rr = hiv.rect.rr + 0.1,
                   hiv.ugc.rr = hiv.ureth.rr + 0.2,
                   hiv.rct.rr = hiv.rect.rr + 0.1,
                   hiv.uct.rr = hiv.ureth.rr + 0.2,
                   hiv.dual.rr = 0.2)

init <- init_msm(nwstats = st,
                 prev.ugc = 0.005, prev.rgc = 0.005,
                 prev.uct = 0.013, prev.rct = 0.013)

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(min.n = njobs, compress = TRUE, outdir = "data/")
