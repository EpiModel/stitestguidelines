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
gctxsuccess <- as.numeric(Sys.getenv("GCTXSUCCESS"))
cttxsuccess <- as.numeric(Sys.getenv("CTTXSUCCESS"))

cat("Array number is ", jobno)
cat("\n fsimno is ", fsimno)

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   tst.rect.sti.rr = 1, #recttest

                   # Correlation
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.5161, #0.513,
                   ugc.tprob = 0.4362, # 0.432
                   rct.tprob = 0.2813, #0.2797,
                   uct.tprob = 0.2195, # 0.2165,
                   syph.tprob = 0, #0.1206,

                   # HIV acquisition
                   hiv.rgc.rr = 1.97, #1.75,
                   hiv.ugc.rr = 1.48, #1.27,
                   hiv.rct.rr = 1.97, #1.75,
                   hiv.uct.rr = 1.48, #1.27,
                   hiv.syph.rr = 1.64,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.82,
                   syph.seco.sympt.prob = 0.90,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.85,
                   syph.seco.sympt.prob.tx = 0.85,
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
                   ept.gc.success = gctxsuccess,
                   ept.ct.success = cttxsuccess,

                   partnercut = 1,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0, #0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0, #0.67,
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
netsim_hpc("est/stimod.burnin.ept.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 truncate.at = 5200,
                 vars =
                   c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                     "ir100.ct", "incid.ct",
                     "incid.sti",
                     "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                     "ir100.sti", "ir100.sti.prep", "ir100.gcct",
                     "incid.gc.hivneg", "incid.gc.hivpos",
                     "incid.ct.hivneg", "incid.ct.hivpos",
                     "ir100.gc.hivneg", "ir100.gc.hivpos",
                     "ir100.ct.hivneg", "ir100.ct.hivpos",
                     "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                     'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                     'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                     "i.prev", "prev.gc", "prev.rgc", "prev.ugc",
                     "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                     "txGC", "txCT",
                     "txGC_asympt", "txCT_asympt", "txSTI", "txSTI_asympt",
                     "recentpartners", "recentpartners.prop",
                     "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
                     "eptTx", "propindexeptElig",
                     "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
                     "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
                     "eptgcctinfectundiaghiv",
                     "gc.timesInf", "ct.timesInf", "sti.timesInf",
                     "eptgcinfecthiv", "eptctinfecthiv",
                     "eptgcctinfecthiv",
                     "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
                     "eptgcctinfecthiv_inst",
                     "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
                     "eptgcctinfectundiaghiv_inst",
                     "eptindexprovided_gc", "eptindexprovided_ct",
                     "eptpartprovided_gc", "eptpartprovided_ct",
                     "eptpartprovided_main", "eptpartprovided_pers",
                     "eptpartprovided_inst", "eptpartuptake_main",
                     "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
                     "eptpartuptake_pers", "eptpartuptake_inst",
                     "eptpartuptake_gc", "eptpartuptake_ct"))
