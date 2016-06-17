
## Test Script for stiPrEP Project

rm(list=ls())
suppressMessages(library(EpiModelHIVmsm))
sourceDir("source/", TRUE)

load("est/nwstats.rda")

param <- param.msm(nwstats = st,
                   testing.pattern = "memoryless",
                   ai.scale = 1.1,
                   riskh.start = 5000,
                   prep.start = 5000,
                   prep.elig.model = "cdc3",
                   prep.class.prob = reallocate_pcp(reall = 0),
                   prep.class.hr = c(1, 0.69, 0.19, 0.05),
                   prep.coverage = 0.4,
                   prep.cov.method = "curr",
                   prep.cov.rate = 1,
                   prep.tst.int = 90,
                   prep.risk.int = 182,

                   rcomp.prob = 0,
                   rcomp.adh.groups = 0:4,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE,

                   rgc.tprob = 0.60,
                   ugc.tprob = 0.48,
                   rct.tprob = 0.40,
                   uct.tprob = 0.32,

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

                   gc.prob.cease = 0,
                   ct.prob.cease = 0,

                   gc.prob.tx = 0.90,
                   ct.prob.tx = 0.85,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,

                   sti.cond.rr = 0.3,

                   hiv.gc.rr = 2,
                   hiv.ct.rr = 2)

init <- init.msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.1,
                 prev.rgc = 0.1,
                 prev.uct = 0.1,
                 prev.rct = 0.1)

control <- control.msm(simno = 1,
                       nsteps = 500,
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
                       stitx.FUN = sti_tx,
                       verbose.FUN = verbose.sti,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN", "test.FUN", "tx.FUN",
                                        "prep.FUN", "progress.FUN", "vl.FUN", "edgescorr.FUN",
                                        "resim_nets.FUN", "disclose.FUN", "acts.FUN", "condoms.FUN",
                                        "riskhist.FUN", "position.FUN", "trans.FUN",
                                        "stitrans.FUN", "stirecov.FUN", "stitx.FUN",
                                        "prev.FUN"))


load("est/fit.rda")
sim <- netsim(est, param, init, control)

plot(sim, y = c("prev.rgc", "prev.ugc", "prev.rct", "prev.uct"), mean.col = 1:4, leg = TRUE)

dat <- initialize.sti(est, param, init, control, s = 1)
for (at in 2:dat$control$nsteps) {
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
  dat <- sti_recov(dat, at)
  dat <- sti_tx(dat, at)
  dat <- prevalence.msm(dat, at)
  cat(at, ".", sep = "")
}

