library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
load("est/nwstats.rda")
anncov <- 0.1
hrcov <- 0.1
anncov <- 0.0
hrcov <- 0.0
annint <- 364
hrint <- 182
partnercutoff <- 1
stiasymptx <- 1
param <- param_msm(nwstats = st,
                   ai.scale = 1.03,

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
                   rgc.tprob = 0.447,
                   ugc.tprob = 0.337,
                   rct.tprob = 0.2025,
                   uct.tprob = 0.1825,
                   syph.tprob = 0.1424,

                   # HIV acquisition
                   hiv.rgc.rr = 1.80292790,
                   hiv.ugc.rr = 1.1989083,
                   hiv.rct.rr = 1.80292790,
                   hiv.uct.rr = 1.1989083,
                   hiv.syph.rr = 1.500918,

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

                   hivdx.syph.sympt.tx.rr = 1.5,

                   partnercut = partnercutoff,
                   stianntest.coverage = anncov,
                   stihighrisktest.coverage = hrcov,
                   prep.coverage = 0,
                   ept.coverage = 0,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,

                   stitest.active.int = annint,
                   sti.highrisktest.int = hrint) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(start = 5201,
                       nsteps = 5720,
                       nsims = 1,
                       ncores = 1,
                       initialize.FUN = reinit_msm,
                       verbose = TRUE)

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = TRUE, verbose = TRUE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE)
