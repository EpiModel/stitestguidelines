
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
library("EpiModel")

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
anngcnegcov <- as.numeric(Sys.getenv("ANNGCNEGCOV"))
anngcposcov <- as.numeric(Sys.getenv("ANNGCPOSCOV"))
annctnegcov <- as.numeric(Sys.getenv("ANNCTNEGCOV"))
annctposcov <- as.numeric(Sys.getenv("ANNCTPOSCOV"))
hrgcnegcov <- as.numeric(Sys.getenv("HRGCNEGCOV"))
hrgcposcov <- as.numeric(Sys.getenv("HRGCPOSCOV"))
hrctnegcov <- as.numeric(Sys.getenv("HRCTNEGCOV"))
hrctposcov <- as.numeric(Sys.getenv("HRCTPOSCOV"))
annint <- as.numeric(Sys.getenv("ANNINT"))
hrint <- as.numeric(Sys.getenv("HRINT"))
partnercutoff <- as.numeric(Sys.getenv("PART"))
stiasymptx <- as.numeric(Sys.getenv("STIASYMPTX"))

cat("Array number is ", jobno)
cat("\n fsimno is ", fsimno)

## Parameters
load("est/nwstats.20k.rda")

param <- param_msm(nwstats = st,

                   # AI Scale
                   ai.scale = 1.061338,
                   ai.scale.pospos = 1.061338,

                   # Probability of rectal testing
                   tst.rect.sti.rr = 1,

                   # Correlation
                   sti.correlation.time = 0,

                   # STI acquisition
                   rgc.tprob = 0.5425, #0.5364416,
                   ugc.tprob = 0.4405, #0.434692,
                   rct.tprob = 0.2480, #0.2493814,
                   uct.tprob = 0.1940, #0.1944415,

                   # HIV acquisition
                   hiv.rgc.rr = 2.175918,
                   hiv.ugc.rr = 1.564797,
                   hiv.rct.rr = 2.175918,
                   hiv.uct.rr = 1.564797,

                   # STI symptom probability
                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.80,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.48,

                   # Asymptomatic probability of treatment
                   gc.asympt.prob.tx = stiasymptx,
                   ct.asympt.prob.tx = stiasymptx,

                   # Partner cutoff
                   partnercut = partnercutoff,

                   # Testing coverage
                   stianntest.gc.hivneg.coverage = anngcnegcov,
                   stianntest.ct.hivneg.coverage = annctnegcov,
                   stianntest.syph.hivneg.coverage = 0,
                   stihighrisktest.gc.hivneg.coverage = hrgcnegcov,
                   stihighrisktest.ct.hivneg.coverage = hrctnegcov,
                   stihighrisktest.syph.hivneg.coverage = 0, #hrsyphnegcov,
                   stianntest.gc.hivpos.coverage = anngcposcov,
                   stianntest.ct.hivpos.coverage = annctposcov,
                   stianntest.syph.hivpos.coverage = 0,
                   stihighrisktest.gc.hivpos.coverage = hrgcposcov,
                   stihighrisktest.ct.hivpos.coverage = hrctposcov,
                   stihighrisktest.syph.hivpos.coverage = 0, #hrsyphposcov,

                   # Other intervention coverage
                   prep.coverage = 0,
                   ept.coverage = 0,

                   # Interventions start
                   prep.start = 7000,
                   stitest.start = 1,
                   ept.start = 7000,

                   # Testing interval
                   stitest.active.int = annint,
                   sti.highrisktest.int = hrint) # adjustable for 3 or 6 months

init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 3120,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       prev.FUN = prevalence_msm_tnt,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/sti.tnt.burnin.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 truncate.at = 2600,
                 vars <- c(
                   # HIV
                   "incid", "hivtests.nprep", "ir100", "i.prev",

                   # GC
                   "incid.rgc", "incid.ugc",
                   "incid.rgc.tttraj1", "incid.ugc.tttraj1",
                   "incid.rgc.tttraj2", "incid.ugc.tttraj2",
                   "ir100.gc", "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                   "prev.gc", "prev.gc.tttraj1", "prev.gc.tttraj2",
                   "GCasympttests", "GCsympttests",
                   "GCasympttests.tttraj1", "GCasympttests.tttraj2",
                   "GCsympttests.tttraj1", "GCsympttests.tttraj2",
                   "txGC", "txGC.tttraj1", "txGC.tttraj2",
                   "txGC_asympt",
                   "tt.traj.gc1", "tt.traj.gc2",

                   # CT
                   "incid.rct", "incid.uct",
                   "incid.rct.tttraj1", "incid.uct.tttraj1",
                   "incid.rct.tttraj2", "incid.uct.tttraj2",
                   "ir100.ct", "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                   "prev.ct", "prev.ct.tttraj1", "prev.ct.tttraj2",
                   "CTasympttests", "CTsympttests",
                   "CTasympttests.tttraj1", "CTasympttests.tttraj2",
                   "CTsympttests.tttraj1", "CTsympttests.tttraj2",
                   "txCT", "txCT.tttraj1", "txCT.tttraj2",
                   "txCT_asympt",
                   "tt.traj.ct1", "tt.traj.ct2",

                   # Combined
                   "incid.gcct", "incid.gcct.tttraj1", "incid.gcct.tttraj2",
                   "ir100.gcct","ir100.ct.tttraj1", "ir100.ct.tttraj2",
                   "prev.STI", "prev.sti.tttraj1", "prev.sti.tttraj2",
                   "stiasympttests", "stisympttests",
                   "stiasympttests.tttraj1", "stiasympttests.tttraj2",
                   "stisympttests.tttraj1", "stisympttests.tttraj2",
                   "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
                   "txSTI_asympt",
                   "tt.traj.sti1", "tt.traj.sti2",

                   # Other
                   "num",
                   'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                   'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg'))
