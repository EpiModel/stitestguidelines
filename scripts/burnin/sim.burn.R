
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

                   # Correlation
                   sti.stitx.correlation = "false",
                   sti.hivdx.correlation = "false",
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.478,
                   ugc.tprob = 0.383,
                   rct.tprob = 0.257,
                   uct.tprob = 0.211,
                   syph.tprob = 0.2525,

                   # HIV acquisition
                   hiv.rgc.rr = 1.90,
                   hiv.ugc.rr = 1.40,
                   hiv.rct.rr = 1.90,
                   hiv.uct.rr = 1.40,
                   hiv.syph.rr = 1.80,

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

                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.43,
                   stianntest.syph.hivneg.coverage = 0.43,
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.60,
                   stianntest.syph.hivpos.coverage = 0.65,
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,
                   # 2014 (Hoots): NG: 46.2% HIV-MSM, 64.1%  HIV+ MSM
                   # 2014 (Hoots): CT: 45.8% HIV-MSM, 62.8%  HIV+ MSM
                   # NHBS syphilis testing (2014 self-report data Qian An): 45% HIV- MSM, 68% HIV+ MSM

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st)

control <- control_msm(simno = fsimno,
                       nsteps = 5200,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = TRUE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
