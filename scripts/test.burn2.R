
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

                   # STI acquisition
                   rgc.tprob = 0.4773,
                   ugc.tprob = 0.3819,
                   rct.tprob = 0.2564,
                   uct.tprob = 0.2091,
                   syph.tprob = 0.2526,

                   # HIV acquisition
                   hiv.rgc.rr = 1.75,
                   hiv.ugc.rr = 1.26,
                   hiv.rct.rr = 1.75,
                   hiv.uct.rr = 1.26,
                   hiv.syph.rr = 1.63,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.70,
                   syph.seco.sympt.prob = 0.85,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.80,
                   syph.seco.sympt.prob.tx = 0.80,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   ept.coverage = 0.5,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.1,
                   stihighrisktest.ct.hivneg.coverage = 0.1,
                   stihighrisktest.syph.hivneg.coverage = 0.1,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.67,
                   stihighrisktest.gc.hivpos.coverage = 0.1,
                   stihighrisktest.ct.hivpos.coverage = 0.1,
                   stihighrisktest.syph.hivpos.coverage = 0.1,

                   prep.start = 7000,
                   stitest.start = 7000,
                   ept.start = 5201,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60,
                   partlist.start = 1)

init <- init_msm(nwstats = st,
                 prev.B = 0.10,
                 prev.W = 0.10,
                 prev.ugc = 0.010,
                 prev.rgc = 0.010,
                 prev.uct = 0.010,
                 prev.rct = 0.010,
                 prev.syph.B = 0.010,
                 prev.syph.W = 0.010)

# init <- init_msm(nwstats = st,
#                  prev.B = 0,
#                  prev.W = 0,
#                  prev.ugc = 0,
#                  prev.rgc = 0,
#                  prev.uct = 0,
#                  prev.rct = 0,
#                  prev.syph.B = 0.05,
#                  prev.syph.W = 0.05)

control <- control_msm(simno = 1,
                       nsteps = 260,
                       nsims = 1,
                       ncores = 1,
                       verbose = TRUE)

# data(est)
sim <- netsim(est, param, init, control)

plot(sim, y = "prev.syph", mean.smooth = FALSE, mean.lwd = 1)
plot(sim, y = "incid.syph", mean.smooth = FALSE, mean.lwd = 1)

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
