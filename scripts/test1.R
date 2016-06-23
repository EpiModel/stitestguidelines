
## Test Script for stiPrEP Project

rm(list=ls())
suppressMessages(library(EpiModelHIVmsm))
EpiModelHPC::sourceDir("source/", TRUE)

devtools::load_all("~/Dropbox/Dev/EpiModelHIVmsm/EpiModelHIVmsm")

# Main Test Script ----------------------------------------------------

load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   ai.scale = 1,

                   riskh.start = 5000,
                   prep.start = 5000,
                   prep.coverage = 0,

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

init <- init_msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.1,
                 prev.rgc = 0.1,
                 prev.uct = 0.1,
                 prev.rct = 0.1)

control <- control_msm(simno = 1,
                       nsteps = 100,
                       nsims = 1,
                       ncores = 1,
                       save.int = 5000,
                       acts.FUN = acts_sti,
                       condoms.FUN = condoms_sti,
                       initialize.FUN = initialize_sti,
                       prep.FUN = prep_sti,
                       prev.FUN = prevalence_sti,
                       riskhist.FUN = riskhist_sti,
                       position.FUN = position_sti,
                       trans.FUN = trans_sti,
                       stitrans.FUN = sti_trans,
                       stirecov.FUN = sti_recov,
                       stitx.FUN = sti_tx,
                       verbose.FUN = verbose_sti,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN", "test.FUN", "tx.FUN",
                                        "prep.FUN", "progress.FUN", "vl.FUN", "edgescorr.FUN",
                                        "resim_nets.FUN", "disclose.FUN", "acts.FUN", "condoms.FUN",
                                        "riskhist.FUN", "position.FUN", "trans.FUN",
                                        "stitrans.FUN", "stirecov.FUN", "stitx.FUN",
                                        "prev.FUN"))


load("est/fit.rda")
sim <- netsim(est, param, init, control)

plot(sim, y = c("prev.rgc", "prev.ugc", "prev.rct", "prev.uct"), mean.col = 1:4, leg = TRUE)



# Testing/Timing ------------------------------------------------------

dat <- initialize_sti(est, param, init, control, s = 1)
for (at in 2:dat$control$nsteps) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- test_sti(dat, at)
  dat <- tx_msm(dat, at)
  dat <- prep_sti(dat, at)
  dat <- progress_msm(dat, at)
  dat <- update_vl_msm(dat, at)
  dat <- edges_correct_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- disclose_msm(dat, at)
  dat <- acts_sti(dat, at)
  dat <- condoms_sti(dat, at)
  dat <- riskhist_sti(dat, at)
  dat <- position_sti(dat, at)
  dat <- trans_sti(dat, at)
  dat <- sti_trans(dat, at)
  dat <- sti_recov(dat, at)
  dat <- sti_tx(dat, at)
  dat <- prevalence_msm(dat, at)
  cat(at, ".", sep = "")
}

library(microbenchmark)
res <- microbenchmark(simnet_msm(dat, at = 2))
summary(res)
