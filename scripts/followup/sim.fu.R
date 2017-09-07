
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
anncov <- as.numeric(Sys.getenv("ANNCOV"))
hrcov <- as.numeric(Sys.getenv("HRCOV"))
annint <- as.numeric(Sys.getenv("ANNINT"))
hrint <- as.numeric(Sys.getenv("HRINT"))
partnercutoff <- as.numeric(Sys.getenv("PART"))
stiasymptx <- as.numeric(Sys.getenv("STIASYMPTX"))

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   ai.scale = 1.03,

                   # Correlation
                   sti.stitx.correlation = "false",
                   sti.hivdx.correlation = "false",
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.4456,
                   ugc.tprob = 0.3341,
                   rct.tprob = 0.1985,
                   uct.tprob = 0.1787,
                   syph.tprob = 0.1464,

                   # HIV acquisition
                   hiv.rgc.rr = 1.80292790,
                   hiv.ugc.rr = 1.1989083,
                   hiv.rct.rr = 1.80292790,
                   hiv.uct.rr = 1.1989083,
                   hiv.syph.rr = 1.62,

                   # HIV transmission
                   hiv.trans.gc.rr = 1.0,
                   hiv.trans.ct.rr = 1.0,
                   hiv.trans.syph.rr = 1.0,

                   syph.prim.sympt.prob.tx = 0.60,
                   syph.seco.sympt.prob.tx = 0.688235,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   syph.prim.asympt.prob.tx = stiasymptx,
                   syph.seco.asympt.prob.tx = stiasymptx,
                   syph.earlat.asympt.prob.tx = stiasymptx,
                   syph.latelat.asympt.prob.tx = stiasymptx,
                   syph.tert.asympt.prob.tx = stiasymptx,
                   gc.asympt.prob.tx = stiasymptx,
                   ct.asympt.prob.tx = stiasymptx,

                   partnercut = partnercutoff,

                   stianntest.gc.hivneg.coverage = anncov, #0.1  + (anncov * 0.1),
                   stianntest.ct.hivneg.coverage = anncov, #0.1 + (anncov * 0.1),
                   stianntest.syph.hivneg.coverage = anncov, #0.1 + (anncov * 0.1),
                   stihighrisktest.gc.hivneg.coverage = hrcov
                   stihighrisktest.ct.hivneg.coverage = hrcov,
                   stihighrisktest.syph.hivneg.coverage = hrcov
                   stianntest.gc.hivpos.coverage = anncov, #0.1 + (anncov * 0.1),
                   stianntest.ct.hivpos.coverage = anncov, #0.1 + (anncov * 0.1),
                   stianntest.syph.hivpos.coverage = anncov, #0.1 + 0.1 + (anncov * 0.1),
                   stihighrisktest.gc.hivpos.coverage = hrcov,
                   stihighrisktest.ct.hivpos.coverage = hrcov,
                   stihighrisktest.syph.hivpos.coverage = hrcov,
                   prep.coverage = 0,
                   ept.coverage = 0,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,

                   stitest.active.int = annint,
                   sti.highrisktest.int = hrint) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 5201,
                       nsteps = 5720,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE)
