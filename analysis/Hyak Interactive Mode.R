library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
load("est/nwstats.rda")
# anncov <- 0.1
# hrcov <- 0.1
# anncov <- 0.0
# hrcov <- 0.0
# annint <- 364
# hrint <- 182
# partnercutoff <- 1
# stiasymptx <- 1

eptcov <- 0.1
prov.main.ong <- 0.5
prov.pers.ong <- 0.4
prov.main.end <- 0.4
prov.pers.end <- 0.3
prov.inst <- 0.2
uptake.main <- 0.8
uptake.pers <- 0.8
uptake.inst <- 0.8
eptint <- 60

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

                   # syph.prim.asympt.prob.tx = stiasymptx,
                   # syph.seco.asympt.prob.tx = stiasymptx,
                   # syph.earlat.asympt.prob.tx = stiasymptx,
                   # syph.latelat.asympt.prob.tx = stiasymptx,
                   # syph.tert.asympt.prob.tx = stiasymptx,
                   # gc.asympt.prob.tx = stiasymptx,
                   # ct.asympt.prob.tx = stiasymptx,

                   ept.provision.partner.main.ong = prov.main.ong,
                   ept.provision.partner.pers.ong = prov.pers.ong,
                   ept.provision.partner.main.end = prov.main.end,
                   ept.provision.partner.pers.end = prov.pers.end,
                   ept.provision.partner.inst = prov.inst,
                   ept.uptake.partner.main = uptake.main,
                   ept.uptake.partner.pers = uptake.pers,
                   ept.uptake.partner.inst = uptake.inst,

                   hivdx.syph.sympt.tx.rr = 1.5,

                   # partnercut = partnercutoff,
                   # stianntest.coverage = anncov,
                   # stihighrisktest.coverage = hrcov,
                   partnercut = 1,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.0,
                   prep.coverage = 0,
                   ept.coverage = eptcov,
                   ept.risk.int = eptint,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(start = 5201,
                       nsteps = 5210,
                       nsims = 3,
                       ncores = 1,
                       initialize.FUN = reinit_msm,
                       verbose = TRUE)

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = TRUE, verbose = TRUE)

process_simfiles(simno = simno, min.n = 3,
                 outdir = "data/", compress = TRUE,
                 delete.sub = TRUE,
                 #truncate.at = 5200,
                 vars =
                   c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                     "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                     "ir100.sti",
                     "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                     "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                     "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                     "ir100.syph.tttraj1", "ir100.syph.tttraj2"))
