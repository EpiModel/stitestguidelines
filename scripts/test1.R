
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

                   rgc.tprob = 0.25,
                   ugc.tprob = 0.25,
                   rct.tprob = 0.25,
                   uct.tprob = 0.25,

                   rgc.sympt.prob = 0.16,
                   ugc.sympt.prob = 0.90,
                   rct.sympt.prob = 0.14,
                   uct.sympt.prob = 0.58,

                   rgc.dur.asympt = 300 / 7,
                   ugc.dur.asympt = 240 / 7,
                   gc.dur.tx = 13 / 7,
                   gc.dur.ntx = 185 / 7,

                   rct.dur.asympt = 497 / 7,
                   uct.dur.asympt = 240 / 7,
                   ct.dur.tx = 14 / 7,
                   ct.dur.ntx = 180 / 7,

                   sti.cond.rr = 0.3,

                   hiv.gc.rr = 2,
                   hiv.ct.rr = 2)

init <- init.msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.05,
                 prev.rgc = 0.05,
                 prev.uct = 0.05,
                 prev.rct = 0.05)

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
                       stirecov.FUN = sti_recov,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN", "test.FUN", "tx.FUN",
                                        "prep.FUN", "progress.FUN", "vl.FUN", "edgescorr.FUN",
                                        "resim_nets.FUN", "disclose.FUN", "acts.FUN", "condoms.FUN",
                                        "riskhist.FUN", "position.FUN", "trans.FUN",
                                        "stitrans.FUN", "stirecov.FUN",
                                        "prev.FUN"))


load("est/fit.10k.rda")
sim <- netsim(est, param, init, control)

sim$epi$prev.rgc
sim$epi$prev.ugc
sim$epi$prev.rct
sim$epi$prev.uct

dat <- initialize.sti(est, param, init, control, s = 1)
at = 2

dat <- aging.msm(dat, at)
dat <- deaths.msm(dat, at)
dat <- births.msm(dat, at)
dat <- test.sti(dat, at)
dat <- tx.msm(dat, at)
dat <- prep.sti(dat, at)
dat <- progress.msm(dat, at)
dat <- update_vl.msm(dat, at)
dat <- edges_correct.msm(dat, at)
dat <- simnet.msm(dat, at)
dat <- disclose.msm(dat, at)
dat <- acts.sti(dat, at)
dat <- condoms.sti(dat, at)
dat <- riskhist.sti(dat, at)
dat <- position.sti(dat, at)
dat <- trans.sti(dat, at)
dat <- sti_trans(dat, at)
# dat <- prevalence.msm(dat, at)
