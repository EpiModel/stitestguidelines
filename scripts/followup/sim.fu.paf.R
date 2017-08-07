
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
                   rgc.tprob = 0.4434,
                   ugc.tprob = 0.3343,
                   rct.tprob = 0.2008,
                   uct.tprob = 0.1790,
                   syph.tprob = 0.1424,

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

                   syph.prim.sympt.prob.tx = 0.60,
                   syph.seco.sympt.prob.tx = 0.688235,
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

                   hivdx.syph.sympt.tx.rr = 1.5,

                   partnercut = 1,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.0,
                   prep.coverage = 0,
                   ept.coverage = 0,

                   prep.start = 7000,
                   stitest.start = 5201,
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
