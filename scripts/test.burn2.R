
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
#devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

# Main Test Script ----------------------------------------------------

# data(st)
load("est/fit.rda")
load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   tst.rect.sti.rr = 1,

                   # STI acquisition
                   rgc.tprob = 0, # 0.4773,
                   ugc.tprob = 0, # 0.3819,
                   rct.tprob = 0, # 0.2564,
                   uct.tprob = 0, # 0.2091,
                   syph.tprob = 0.06,

                   # HIV acquisition
                   hiv.rgc.rr = 1.75,
                   hiv.ugc.rr = 1.26,
                   hiv.rct.rr = 1.75,
                   hiv.uct.rr = 1.26,
                   hiv.syph.rr = 1.63,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.75,
                   syph.seco.sympt.prob = 0.90,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.85,
                   syph.seco.sympt.prob.tx = 0.85,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   ept.coverage = 0.5,
                   stianntest.gc.hivneg.coverage = 0.40, #0.44,
                   stianntest.ct.hivneg.coverage = 0.40, #0.44,
                   stianntest.syph.hivneg.coverage = 0.40, #0.45
                   stihighrisktest.gc.hivneg.coverage = 0.1,
                   stihighrisktest.ct.hivneg.coverage = 0.1,
                   stihighrisktest.syph.hivneg.coverage = 0.1,
                   stianntest.gc.hivpos.coverage = 0.55, #0.61,
                   stianntest.ct.hivpos.coverage = 0.55, #0.61,
                   stianntest.syph.hivpos.coverage = 0.60, #0.67
                   stihighrisktest.gc.hivpos.coverage = 0.1,
                   stihighrisktest.ct.hivpos.coverage = 0.1,
                   stihighrisktest.syph.hivpos.coverage = 0.1,

                   # Condoms
                   cond.main.BB.prob = 0.21, # 0.21,
                   cond.main.BW.prob = 0.21, # 0.21,
                   cond.main.WW.prob = 0.21, # 0.21,
                   cond.pers.always.prob = 0.216, # 0.216,
                   cond.pers.BB.prob = 0.26, # 0.26,
                   cond.pers.BW.prob = 0.26, # 0.26,
                   cond.pers.WW.prob = 0.26, # 0.26,
                   cond.inst.always.prob = 0.326, # 0.326,
                   cond.inst.BB.prob = 0.27, # 0.27,
                   cond.inst.BW.prob = 0.27, # 0.27,
                   cond.inst.WW.prob = 0.27, # 0.27,

                   prep.start = 7000,
                   stitest.start = 2601,
                   ept.start = 5201,

                   #partlist.start = 1,
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60
)

init <- init_msm(nwstats = st,
                 prev.B = 0.10,
                 prev.W = 0.10,
                 prev.ugc = 0.010,
                 prev.rgc = 0.010,
                 prev.uct = 0.010,
                 prev.rct = 0.010,
                 prev.syph.B = 0.010,
                 prev.syph.W = 0.010)

control <- control_msm(simno = 1,
                       nsteps = 2600,
                       nsims = 1,
                       ncores = 1,
                       verbose = TRUE)

# data(est)
sim <- netsim(est, param, init, control)

plot(sim, y = "prev.syph", mean.smooth = FALSE, mean.lwd = 1)
plot(sim, y = "prev.primsecosyph", mean.smooth = FALSE, mean.lwd = 1)
plot(sim, y = "incid.syph", mean.smooth = FALSE, mean.lwd = 1)
plot(sim, y = "ir100.syph", mean.smooth = FALSE, mean.lwd = 1)
plot(sim, y = "i.prev", mean.smooth = FALSE, mean.lwd = 1)

# Empirical values
# 2014 (Hoots): NG: 46.2% HIV-MSM, 64.1%  HIV+ MSM
# 2014 (Hoots): CT: 45.8% HIV-MSM, 62.8%  HIV+ MSM
# NHBS syphilis testing (2014 self-report data Qian An):
# 45% HIV- MSM, 68% HIV+ MSM

par(mfrow = c(2, 2))
plot(sim, y = "test.ct.12mo.hivdiag")
title("CT + HIV diag test in past 12 mos")
abline(h = 0.628)
plot(sim, y = "test.gc.12mo.hivdiag")
abline(h = 0.641)
title("NG + HIV diag test in past 12 mos")
plot(sim, y = "test.syph.12mo.hivdiag")
abline(h = 0.68)
title("Syph + HIV diag test in past 12 mos")
abline(h = 0.8)

par(mfrow = c(2, 2))
plot(sim, y = "test.ct.12mo.nonhivdiag")
abline(h = 0.458)
title("CT + non-HIV diag test in past 12 mos")
plot(sim, y = "test.gc.12mo.nonhivdiag")
abline(h = 0.462)
title("NG + non-HIV diag test in past 12 mos")
plot(sim, y = "test.syph.12mo.nonhivdiag")
abline(h = 0.45)
title("Syph + non-HIV diag test in past 12 mos")

par(mfrow = c(2, 2))
plot(sim, y = "ir100")
title("HIV Incidence/ 100 PY")
plot(sim, y = "ir100.gc")
abline(h = 3.5)
title("GC Incidence/ 100 PY")
plot(sim, y = "ir100.ct")
abline(h = 5.6)
title("CT Incidence/ 100 PY")
plot(sim, y = "ir100.syph")
abline(h = 2.6)
title("Syph Incidence/ 100 PY")

par(mfrow = c(2, 2))
plot(sim, y = "i.prev")
title("HIV Prevalence")
abline(h = 0.15)
plot(sim, y = "prev.syph")
abline(h = 0.012)
title("Syph all-stages prevalence")
plot(sim, y = "prev.primsecosyph")
abline(h = 0.003)
title("P&S Syph prevalence")


df <- as.data.frame(sim)
names(df)
mean(tail(df$prev.syph, 10)) # 0.05248137
mean(tail(df$incid.syph, 25)) # 2.28

# Testing/Timing ------------------------------------------------------

# control$bi.mods
debug(syph_progress_msm)

dat <- initialize_msm(est, param, init, control, s = 1)

for (at in 2:50) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- hiv_test_msm(dat, at)
  dat <- sti_test_msm(dat, at)
  dat <- hiv_tx_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- hiv_progress_msm(dat, at)
  dat <- syph_progress_msm(dat, at)
  dat <- hiv_vl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- hiv_disclose_msm(dat, at)
  dat <- part_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- riskhist_prep_msm(dat, at)
  dat <- riskhist_stitest_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- hiv_trans_msm(dat, at)
  dat <- sti_trans_msm(dat, at)
  dat <- sti_recov_msm(dat, at)
  dat <- sti_tx_msm(dat, at)
  dat <- sti_ept_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  cat(at, " - ", sep = "")
}

str(dat$temp$part.list)
head(dat$temp$part.list)
pl <- dat$temp$part.list
# pl <- pl[pl[, "ptype"] == 1, ]
table(pl[, 6])

ids <- c(pl[, 1], pl[, 2])
tabs <- tabulate(ids)
summary(tabs)
hist(tabs)
mean(tabs > 10)
