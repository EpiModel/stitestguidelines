
library("methods")
suppressMessages(library("EpiModelHIV"))

load("est/nwstats.rda")

load("scripts/burnin/abc/abc.avg.parms.1pct.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,

                   rgc.tprob = rgc.tprob,
                   ugc.tprob = ugc.tprob,
                   rct.tprob = rct.tprob,
                   uct.tprob = uct.tprob,

                   rgc.sympt.prob = rgc.sympt.prob,
                   ugc.sympt.prob = ugc.sympt.prob,
                   rct.sympt.prob = rct.sympt.prob,
                   uct.sympt.prob = uct.sympt.prob,

                   rgc.dur.asympt = rgc.dur.asympt,
                   ugc.dur.asympt = ugc.dur.asympt,

                   rct.dur.asympt = rct.dur.asympt,
                   uct.dur.asympt = uct.dur.asympt,

                   gc.prob.cease = prob.cease,
                   ct.prob.cease = prob.cease,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,
                   prep.continue.stand.tx = TRUE,

                   hiv.rgc.rr = hiv.rect.rr,
                   hiv.ugc.rr = hiv.ureth.rr,
                   hiv.rct.rr = hiv.rect.rr,
                   hiv.uct.rr = hiv.ureth.rr,
                   hiv.dual.rr = 0)

init <- init_msm(st)

control <- control_msm(simno = 1,
                       start = 2601, nsteps = 2700,
                       nsims = 1, ncores = 1,
                       initialize.FUN = reinit_msm)

load("est/stimod.burnin.rda")
sim2 <- netsim(sim, param, init, control)


# Testing/Timing ------------------------------------------------------

debug(sti_tx)

load("est/stimod.burnin.rda")
dat <- reinit_msm(sim, param, init, control, s = 1)

for (at in 2601:2700) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- test_msm(dat, at)
  dat <- tx_msm(dat, at)
  dat <- prep_sti(dat, at)
  dat <- progress_msm(dat, at)
  dat <- vl_msm(dat, at)
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
  cat(at, ".", sep = "")
}

undebug(prep_sti)
debug(riskhist_sti)

