
# system("scp scripts/followup/*.fu.[Rs]* hyak:/gscratch/csde/sjenness/sti")
# system("scp source/*.* hyak:/gscratch/csde/sjenness/sti/source")

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
sourceDir("source/", FALSE)

## Environmental Arguments
# args <- commandArgs(trailingOnly = TRUE)
# simno <- args[1]
# jobno <- args[2]

cov <- 0.5

## Parameters
# fsimno <- paste(simno, jobno, sep = ".")
fsimno = 1
load("est/nwstats.rda")

load("scripts/burnin/abc.parms.v1.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,
                   ai.scale = 1,

                   prep.coverage = cov,
                   prep.start = 2601,

                   rcomp.prob = 0,
                   rcomp.adh.groups = 0:3,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE,

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
                   gc.dur.tx = 2,
                   gc.dur.ntx = NULL,

                   rct.dur.asympt = rct.dur.asympt,
                   uct.dur.asympt = uct.dur.asympt,
                   ct.dur.tx = 2,
                   ct.dur.ntx = NULL,

                   gc.prob.cease = 0,
                   ct.prob.cease = 0,

                   gc.sympt.prob.tx = 0.90,
                   ct.sympt.prob.tx = 0.85,
                   gc.asympt.prob.tx = 0,
                   ct.asympt.prob.tx = 0,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,
                   prep.continue.stand.tx = TRUE,

                   sti.cond.rr = 0.3,

                   hiv.rgc.rr = hiv.rect.rr,
                   hiv.ugc.rr = hiv.ureth.rr,
                   hiv.rct.rr = hiv.rect.rr,
                   hiv.uct.rr = hiv.ureth.rr,
                   hiv.dual.rr = 0)

init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 3120,
                       nsims = 1,
                       ncores = 1,
                       acts.FUN = acts_sti,
                       condoms.FUN = condoms_sti,
                       initialize.FUN = reinit_msm,
                       prep.FUN = prep_sti,
                       prev.FUN = prevalence_sti,
                       riskhist.FUN = riskhist_sti,
                       position.FUN = position_sti,
                       trans.FUN = trans_sti,
                       stitrans.FUN = sti_trans,
                       stirecov.FUN = sti_recov,
                       stitx.FUN = sti_tx,
                       verbose.FUN = verbose_sti,
                       verbose.int = 1,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN",
                                        "test.FUN", "tx.FUN", "prep.FUN",
                                        "progress.FUN", "vl.FUN",
                                        "resim_nets.FUN", "disclose.FUN",
                                        "acts.FUN", "condoms.FUN", "riskhist.FUN",
                                        "position.FUN", "trans.FUN", "stitrans.FUN",
                                        "stirecov.FUN", "stitx.FUN", "prev.FUN"))


load("est/stimod.burnin.rda")
sim2 <- netsim(sim, param, init, control)

debug(sti_tx)


# Testing/Timing ------------------------------------------------------

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

