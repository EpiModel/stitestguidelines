
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- Sys.getenv("SIMNO")
jobno <- Sys.getenv("PBS_ARRAYID")
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

cov <- as.numeric(Sys.getenv("COV"))
prstiint <- as.numeric(Sys.getenv("PSTIINT"))
rc <- as.numeric(Sys.getenv("RC"))
asympt <- as.numeric(Sys.getenv("ASYMPT"))

## Parameters
load("est/nwstats.rda")

load("est/abc.avg.parms.1pct.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,

                   prep.coverage = cov,
                   prep.start = 2601,

                   rcomp.prob = rc,
                   rcomp.adh.groups = 2:3,

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

                   gc.asympt.prob.tx = asympt,
                   ct.asympt.prob.tx = asympt,

                   prep.sti.screen.int = prstiint,

                   hiv.rgc.rr = hiv.rect.rr,
                   hiv.ugc.rr = hiv.ureth.rr,
                   hiv.rct.rr = hiv.rect.rr,
                   hiv.uct.rr = hiv.ureth.rr)

init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 3120,
                       nsims = 50,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/stimod.mean1pct.rda", param, init, control,
           compress = FALSE, verbose = FALSE)

process_simfiles(min.n = njobs, outdir = "data/", compress = "xz")
