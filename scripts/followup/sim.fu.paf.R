
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
syphacq <-as.numeric(Sys.getenv("SYPHACQ"))
rectacq <- as.numeric(Sys.getenv("RECTACQ"))
urethacq <- as.numeric(Sys.getenv("URETHACQ"))
syphtrans <- as.numeric(Sys.getenv("SYPHTRANS"))
gctrans <- as.numeric(Sys.getenv("GCTRANS"))
cttrans <- as.numeric(Sys.getenv("CTTRANS"))
gccttrans <- as.numeric(Sys.getenv("GCCTTRANS"))
ctsyphtrans <- as.numeric(Sys.getenv("CTSYPHTRANS"))
gcsyphtrans <- as.numeric(Sys.getenv("GCSYPHTRANS"))
allstitrans <- as.numeric(Sys.getenv("ALLSTITRANS"))

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,

                   # Correlation
                   sti.stitx.correlation = "false",
                   sti.hivdx.correlation = "false",
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.4773,
                   ugc.tprob = 0.3819,
                   rct.tprob = 0.2564,
                   uct.tprob = 0.2091,
                   syph.tprob = 0.2553,

                   # HIV acquisition
                   hiv.rgc.rr = rectacq,
                   hiv.ugc.rr = urethacq,
                   hiv.rct.rr = rectacq,
                   hiv.uct.rr = urethacq,
                   hiv.syph.rr = syphacq,

                   # HIV transmission
                   hiv.trans.gc.rr = gctrans,
                   hiv.trans.ct.rr = cttrans,
                   hiv.trans.syph.rr = syphtrans,
                   hiv.trans.gc.ct.rr = gccttrans,
                   hiv.trans.gc.syph.rr = gcsyphtrans,
                   hiv.trans.ct.syph.rr = ctsyphtrans,
                   hiv.trans.allsti.rr = allstitrans,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.70,
                   syph.seco.sympt.prob = 0.85,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.80,
                   syph.seco.sympt.prob.tx = 0.80,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   syph.prim.asympt.prob.tx = 1.0,
                   syph.seco.asympt.prob.tx = 1.0,
                   syph.earlat.asympt.prob.tx = 1.0,
                   syph.latelat.asympt.prob.tx = 1.0,
                   syph.tert.asympt.prob.tx = 1.0,
                   gc.asympt.prob.tx = 1.0,
                   ct.asympt.prob.tx = 1.0,

                   partnercut = 1,
                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.67,
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,
                   prep.coverage = 0,

                   prep.start = 7000,
                   stitest.start = 7000,
                   ept.start = 7000,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months


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
