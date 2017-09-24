
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
anngcnegcov <- as.numeric(Sys.getenv("ANNGCNEGCOV"))
anngcposcov <- as.numeric(Sys.getenv("ANNGCPOSCOV"))
annctnegcov <- as.numeric(Sys.getenv("ANNCTNEGCOV"))
annctposcov <- as.numeric(Sys.getenv("ANNCTPOSCOV"))
annsyphnegcov <- as.numeric(Sys.getenv("ANNSYPHNEGCOV"))
annsyphposcov <- as.numeric(Sys.getenv("ANNSYPHPOSCOV"))
hrgcnegcov <- as.numeric(Sys.getenv("HRGCNEGCOV"))
hrgcposcov <- as.numeric(Sys.getenv("HRGCPOSCOV"))
hrctnegcov <- as.numeric(Sys.getenv("HRCTNEGCOV"))
hrctposcov <- as.numeric(Sys.getenv("HRCTPOSCOV"))
hrsyphnegcov <- as.numeric(Sys.getenv("HRSYPHNEGCOV"))
hrsyphposcov <- as.numeric(Sys.getenv("HRSYPHPOSCOV"))
annint <- as.numeric(Sys.getenv("ANNINT"))
hrint <- as.numeric(Sys.getenv("HRINT"))
partnercutoff <- as.numeric(Sys.getenv("PART"))
stiasymptx <- as.numeric(Sys.getenv("STIASYMPTX"))

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
                   hiv.rgc.rr = 1.78,
                   hiv.ugc.rr = 1.29,
                   hiv.rct.rr = 1.78,
                   hiv.uct.rr = 1.29,
                   hiv.syph.rr = 1.66,

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

                   syph.prim.asympt.prob.tx = stiasymptx,
                   syph.seco.asympt.prob.tx = stiasymptx,
                   syph.earlat.asympt.prob.tx = stiasymptx,
                   syph.latelat.asympt.prob.tx = stiasymptx,
                   syph.tert.asympt.prob.tx = stiasymptx,
                   gc.asympt.prob.tx = stiasymptx,
                   ct.asympt.prob.tx = stiasymptx,

                   partnercut = partnercutoff,

                   stianntest.gc.hivneg.coverage = anngcnegcov,
                   stianntest.ct.hivneg.coverage = annctnegcov,
                   stianntest.syph.hivneg.coverage = annsyphnegcov,
                   stihighrisktest.gc.hivneg.coverage = hrgcnegcov,
                   stihighrisktest.ct.hivneg.coverage = hrctnegcov,
                   stihighrisktest.syph.hivneg.coverage = hrsyphnegcov,
                   stianntest.gc.hivpos.coverage = anngcposcov,
                   stianntest.ct.hivpos.coverage = annctposcov,
                   stianntest.syph.hivpos.coverage = annsyphposcov,
                   stihighrisktest.gc.hivpos.coverage = hrgcposcov,
                   stihighrisktest.ct.hivpos.coverage = hrctposcov,
                   stihighrisktest.syph.hivpos.coverage = hrsyphposcov,
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
