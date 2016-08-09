
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
load("abc.avg.parms.1pct.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,

                   rgc.tprob = rgc.tprob,
                   ugc.tprob = ugc.tprob,
                   rct.tprob = rct.tprob,
                   uct.tprob = uct.tprob,

                   rgc.sympt.prob = rgc.sympt.prob,
                   ugc.sympt.prob = ugc.sympt.prob,
                   rct.sympt.prob = rct.sympt.prob,
                   uct.sympt.prob = uct.sympt.prob,

                   rgc.dur.asympt = rgc.dur.asympt,
                   ugc.dur.asympt = ugc.dur.asympt,

                   rct.dur.asympt = rct.dur.asympt,
                   uct.dur.asympt = uct.dur.asympt,

                   gc.prob.cease = prob.cease,
                   ct.prob.cease = prob.cease,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,
                   prep.continue.stand.tx = TRUE,

                   hiv.rgc.rr = hiv.rect.rr,
                   hiv.ugc.rr = hiv.ureth.rr,
                   hiv.rct.rr = hiv.rect.rr,
                   hiv.uct.rr = hiv.ureth.rr,
                   hiv.dual.rr = 0)

init <- init_msm(nwstats = st)

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control, cp.save.int = NULL,
            save.min = TRUE, save.max = TRUE)

process_simfiles(min.n = njobs, outdir = "data/")
