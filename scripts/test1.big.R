
## Test Script for stiPrEP Project

rm(list=ls())
suppressMessages(library(EpiModelHIV))
sourceDir("source/", TRUE)

# devtools::load_all("~/Dropbox/Dev/EpiModelHIVmsm/EpiModelHIVmsm")

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

                   hiv.rgc.rr = 2.2,
                   hiv.ugc.rr = 1.5,
                   hiv.rct.rr = 2.2,
                   hiv.uct.rr = 1.5)

init <- init_msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.1,
                 prev.rgc = 0.1,
                 prev.uct = 0.1,
                 prev.rct = 0.1)

control <- control_msm(simno = 1,
                       nsteps = 5200,
                       nsims = 16,
                       ncores = 16,
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
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN",
                                        "test.FUN", "tx.FUN", "prep.FUN",
                                        "progress.FUN", "vl.FUN",
                                        "resim_nets.FUN", "disclose.FUN",
                                        "acts.FUN", "condoms.FUN", "riskhist.FUN",
                                        "position.FUN", "trans.FUN", "stitrans.FUN",
                                        "stirecov.FUN", "stitx.FUN", "prev.FUN"))


load("est/fit.rda")
sim <- netsim(est, param, init, control)
save(sim, file = "data.rda", compress = "xz")

plot(sim, y = c("prev.rgc", "prev.ugc", "prev.rct", "prev.uct"),
     mean.col = 1:4, leg = TRUE, qnts = 1, qnts.col = 1:4)
#
#
#
# # Testing/Timing ------------------------------------------------------
#
# dat <- initialize_sti(est, param, init, control, s = 1)
#
# for (at in 2:100) {
#   dat <- aging_msm(dat, at)       ## <1 ms
#   dat <- deaths_msm(dat, at)      ## 8 ms
#   dat <- births_msm(dat, at)      ## 9 ms
#   dat <- test_msm(dat, at)        ## 2 ms
#   dat <- tx_msm(dat, at)          ## 3 ms
#   dat <- prep_sti(dat, at)        ## 6 ms
#   dat <- progress_msm(dat, at)    ## 2 ms
#   dat <- vl_msm(dat, at)          ## 3 ms
#   dat <- simnet_msm(dat, at)      ## 92 ms
#   dat <- disclose_msm(dat, at)    ## 1 ms
#   dat <- acts_sti(dat, at)        ## 1 ms
#   dat <- condoms_sti(dat, at)     ## 2 ms
#   dat <- riskhist_sti(dat, at)    ## 66 ms
#   dat <- position_sti(dat, at)    ## 1 ms
#   dat <- trans_sti(dat, at)       ## 1 ms
#   dat <- sti_trans(dat, at)       ## 4 ms
#   dat <- sti_recov(dat, at)       ## 3 ms
#   dat <- sti_tx(dat, at)          ## 2 ms
#   dat <- prevalence_msm(dat, at)  ## 1 ms
#   cat(at, ".", sep = "")
# }
#
# library(microbenchmark)
#
# res <- microbenchmark(f(dat, at = 101))
# summary(res, unit = "ms")

