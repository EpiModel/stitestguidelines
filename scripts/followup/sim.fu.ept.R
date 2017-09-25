## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
eptcov <- as.numeric(Sys.getenv("EPTCOV"))
eptint <- as.numeric(Sys.getenv("EPTINT"))
prov.main.ong <- as.numeric(Sys.getenv("PROVMAINONG"))
prov.pers.ong <- as.numeric(Sys.getenv("PROVPERSONG"))
prov.main.end <- as.numeric(Sys.getenv("PROVMAINEND"))
prov.pers.end <- as.numeric(Sys.getenv("PROVPERSEND"))
prov.inst <- as.numeric(Sys.getenv("PROVINST"))
uptake.main <- as.numeric(Sys.getenv("UPTAKEMAIN"))
uptake.pers <- as.numeric(Sys.getenv("UPTAKEPERS"))
uptake.inst <- as.numeric(Sys.getenv("UPTAKEINST"))

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

                   # HIV transmission
                   hiv.trans.gc.rr = 1.0,
                   hiv.trans.ct.rr = 1.0,
                   hiv.trans.syph.rr = 1.0,

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

                   # EPT
                   ept.provision.partner.main.ong = prov.main.ong,
                   ept.provision.partner.pers.ong = prov.pers.ong,
                   ept.provision.partner.main.end = prov.main.end,
                   ept.provision.partner.pers.end = prov.pers.end,
                   ept.provision.partner.inst = prov.inst,
                   ept.uptake.partner.main = uptake.main,
                   ept.uptake.partner.pers = uptake.pers,
                   ept.uptake.partner.inst = uptake.inst,

                   partnercut = 1,
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

                   ept.coverage = eptcov,
                   ept.risk.int = eptint,

                   prep.start = 7000,
                   stitest.start = 7000,
                   ept.start = 5201,

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
