
## Test Script for stiPrEP Project

library(EpiModelHIVmsm)
sourceDir("source/", TRUE)

load("est/nwstats.10k.rda")

param <- param.msm(nwstats = st,
                   testing.pattern = "interval",
                   ai.scale = 1.323,
                   riskh.start = 1,
                   prep.start = 30,
                   prep.elig.model = "cdc3",
                   prep.class.prob = reallocate_pcp(reall = 0),
                   prep.class.hr = c(1, 0.69, 0.19, 0.05),
                   prep.coverage = 0.5,
                   prep.cov.method = "curr",
                   prep.cov.rate = 1,
                   prep.tst.int = 90,
                   prep.risk.int = 182,

                   rcomp.prob = 0.5,
                   rcomp.adh.groups = 0:4,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE,

                   gc.rt.patp = 0.25,
                   gc.ur.patp = 0.25,
                   ct.rt.patp = 0.25,
                   ct.ur.patp = 0.25,

                   sti.cond.rr = 0.3,

                   hiv.gc.rr = 2,
                   hiv.ct.rr = 2)

init <- init.msm(nwstats = st, prev.B = 0.253, prev.W = 0.253,
                 prev.ur.gc = 0.05,
                 prev.rt.gc = 0.05,
                 prev.ur.ct = 0.05,
                 prev.rt.ct = 0.05)

control <- control.msm(simno = 1,
                       nsteps = 100,
                       nsims = 1,
                       ncores = 1,
                       save.int = 5000,
                       verbose.int = 1,
                       save.other = c("attr", "temp", "riskh", "el", "p"),
                       acts.FUN = acts.sti,
                       condoms.FUN = condoms.sti,
                       initialize.FUN = initialize.sti,
                       prep.FUN = prep.sti,
                       prev.FUN = prevalence.sti,
                       riskhist.FUN = riskhist.sti,
                       position.FUN = position.sti,
                       trans.FUN = trans.sti,
                       test.FUN = test.sti,
                       stitrans.FUN = sti_trans,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN", "test.FUN", "tx.FUN",
                                        "prep.FUN", "progress.FUN", "vl.FUN", "edgescorr.FUN",
                                        "resim_nets.FUN", "disclose.FUN", "acts.FUN", "condoms.FUN",
                                        "riskhist.FUN", "position.FUN", "trans.FUN", "stitrans.FUN", "prev.FUN"))




load("est/fit.10k.rda")
sim <- netsim(est, param, init, control)

