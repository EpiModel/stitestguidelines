
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
prov.main.ong <- as.numeric(Sys.getenv("PROVMAINONG"))
prov.pers.ong <- as.numeric(Sys.getenv("PROVPERSONG"))
prov.main.end <- as.numeric(Sys.getenv("PROVMAINEND"))
prov.pers.end <- as.numeric(Sys.getenv("PROVPERSEND"))
prov.inst <- as.numeric(Sys.getenv("PROVINST"))
uptake.main <- as.numeric(Sys.getenv("UPTAKEMAIN"))
uptake.pers <- as.numeric(Sys.getenv("UPTAKEPERS"))
uptake.inst <- as.numeric(Sys.getenv("UPTAKEINST"))
eptint <- as.numeric(Sys.getenv("EPTINT"))

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

                   # EPT
                   ept.provision.partner.main.ong = prov.main.ong,
                   ept.provision.partner.pers.ong = prov.pers.ong,
                   ept.provision.partner.main.end = prov.main.end,
                   ept.provision.partner.pers.end = prov.pers.end,
                   ept.provision.partner.inst = prov.inst,
                   ept.uptake.partner.main = uptake.main,
                   ept.uptake.partner.pers = uptake.pers,
                   ept.uptake.partner.inst = uptake.inst,

                   hivdx.syph.sympt.tx.rr = 1.5,

                   partnercut = 1,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.0,
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
