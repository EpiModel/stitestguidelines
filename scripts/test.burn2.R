
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

# Main Test Script ----------------------------------------------------

# data(st)
load("est/fit.rda")
load("est/nwstats.rda")

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
                   rgc.tprob = 0.4434,
                   ugc.tprob = 0.3343,
                   rct.tprob = 0.2008,
                   uct.tprob = 0.1790,
                   syph.tprob = 0.2533,

                   # HIV acquisition
                   hiv.rgc.rr = 1.80292790,
                   hiv.ugc.rr = 1.1989083,
                   hiv.rct.rr = 1.80292790,
                   hiv.uct.rr = 1.1989083,
                   hiv.syph.rr = 1.62,

                   # HIV transmission
                   hiv.trans.gc.rr = 1,
                   hiv.trans.ct.rr = 1,
                   hiv.trans.syph.rr = 1,

                   syph.prim.sympt.prob.tx = 0.60,
                   syph.seco.sympt.prob.tx = 0.688235,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   syph.prim.asympt.prob.tx = 1,
                   syph.seco.asympt.prob.tx = 1,
                   syph.earlat.asympt.prob.tx = 1,
                   syph.latelat.asympt.prob.tx = 1,
                   syph.tert.asympt.prob.tx = 1,

                   hivdx.syph.sympt.tx.rr = 1.5,

                   prep.coverage = 0.0,
                   ept.coverage = 0.0,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.1,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,

                   stitest.elig.model = "sti",
                   sti.cond.rr = 0.3,

                   cond.rr.BB = 1,
                   cond.rr.BW = 1,
                   cond.rr.WW = 1,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60,
                   partlist.start = 1)

# init <- init_msm(nwstats = st,
#                  prev.B = 0.10,
#                  prev.W = 0.10,
#                  prev.ugc = 0.010,
#                  prev.rgc = 0.010,
#                  prev.uct = 0.010,
#                  prev.rct = 0.010,
#                  prev.syph.B = 0.010,
#                  prev.syph.W = 0.010)

init <- init_msm(nwstats = st,
                 prev.B = 0,
                 prev.W = 0,
                 prev.ugc = 0,
                 prev.rgc = 0,
                 prev.uct = 0,
                 prev.rct = 0,
                 prev.syph.B = 0.05,
                 prev.syph.W = 0.05)

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
