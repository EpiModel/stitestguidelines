
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- Sys.getenv("SIMNO")
jobno <- Sys.getenv("PBS_ARRAYID")
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

## Parameters
load("est/nwstats.rda")
#load("est/abc.syph.parms.rda")
#for (i in seq_along(mean.p)) {
#    assign(names(mean.p)[i], unname(mean.p[i]))
#}

param <- param_msm(nwstats = st,
                   ai.scale = 1.04,

                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.00015625599,

                   # STI acquisition
                   rgc.tprob = 0.442,
                   ugc.tprob = 0.332,
                   rct.tprob = 0.1995,
                   uct.tprob = 0.1795,
                   rsyph.tprob = 0.1495,
                   usyph.tprob = 0.1295,

                   # HIV acquisition
                   hiv.rgc.rr = 2.00,
                   hiv.ugc.rr = 1.20,
                   hiv.rct.rr = 2.00,
                   hiv.uct.rr = 1.20,
                   hiv.rsyph.rr = 2.00,
                   hiv.usyph.rr = 1.20,

                   # HIV transmission
                   hiv.trans.gc.rr = 1,
                   hiv.trans.ct.rr = 1,
                   hiv.trans.syph.rr = 1,

                   syph.prim.sympt.prob.tx = 0.60,
                   syph.seco.sympt.prob.tx = 0.688235,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   syph.prim.asympt.prob.tx = 1,
                   syph.seco.asympt.prob.tx = 1,
                   syph.earlat.asympt.prob.tx = 1,
                   syph.latelat.asympt.prob.tx = 1,
                   syph.tert.asympt.prob.tx = 1,

                   hivdx.syph.sympt.tx.rr = 1.5,

                   prep.coverage = 0.0,
                   ept.coverage = 0.0,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.1,

                   prep.start = 5000,
                   stitest.start = 2601,
                   ept.start = 5000,

                   stitest.elig.model = "sti",

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.B = 0.13,
                 prev.W = 0.13,
                 prev.ugc = 0.010,
                 prev.rgc = 0.010,
                 prev.uct = 0.015,
                 prev.rct = 0.015,
                 prev.syph.B = 0.015,
                 prev.syph.W = 0.015)

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = TRUE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
