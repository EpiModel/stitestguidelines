
rm(list=ls())
library("methods")
suppressMessages(library("EpiModelHIV"))
# devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

load("est/nwstats.rda")
param <- param_msm(nwstats = st,

                   prep.start = 2601,
                   prep.coverage = 0.4)

init <- init_msm(st)

control <- control_msm(simno = 1,
                       start = 2601, nsteps = 2700,
                       nsims = 1, ncores = 1,
                       initialize.FUN = reinit_msm)

load("est/stimod.burnin.rda")
sim2 <- netsim(sim, param, init, control)

colMeans(sim2$epi$prop.CT.asympt.tx, na.rm = TRUE)
colMeans(sim2$epi$prop.GC.asympt.tx, na.rm = TRUE)
colMeans(sim2$epi$prop.rGC.tx, na.rm = TRUE)
colMeans(sim2$epi$prop.rCT.tx, na.rm = TRUE)


# Testing/Timing ------------------------------------------------------

control$bi.mods

debug(sti_tx)
# debug(sti_recov)
# debug(prevalence_msm)

load("est/stimod.burnin.rda")
dat <- reinit_msm(sim, param, init, control, s = 1)

for (at in 2601:2650) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- test_msm(dat, at)
  dat <- tx_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- progress_msm(dat, at)
  dat <- vl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- disclose_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- riskhist_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- trans_msm(dat, at)
  dat <- sti_trans(dat, at)
  dat <- sti_recov(dat, at)
  dat <- sti_tx(dat, at)
  dat <- prevalence_msm(dat, at)
  cat(at, ".", sep = "")
}

undebug(prep_msm)
debug(sti_tx)
