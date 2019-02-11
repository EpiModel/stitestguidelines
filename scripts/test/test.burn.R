rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

# Main Test Script ----------------------------------------------------

data(st)
data(est)
param <- param_msm(nwstats = st,
                   ai.scale = 1.04,

                   # STI acquisition
                   rgc.tprob = 0.4773,
                   ugc.tprob = 0.3819,
                   rct.tprob = 0.2564,
                   uct.tprob = 0.2091,
                   syph.tprob = 0.2533,

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

                   prep.start = 7000,
                   stitest.start = 5201,

                   stitest.elig.model = "all",

                   prep.coverage = 0.0,
                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.10,
                   stihighrisktest.ct.hivneg.coverage = 0.10,
                   stihighrisktest.syph.hivneg.coverage = 0.10,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.67,
                   stihighrisktest.gc.hivpos.coverage = 0.10,
                   stihighrisktest.ct.hivpos.coverage = 0.10,
                   stihighrisktest.syph.hivpos.coverage = 0.10,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months
init <- init_msm(nwstats = st)

control <- control_msm(simno = 1,
                       nsteps = 100,
                       nsims = 1,
                       ncores = 1)


sim <- netsim(est, param, init, control)

df <- as.data.frame(sim)
names(df)
df$incid.gc
df$ir100.gc
df$prev.gc.dual
df$prev.ct.dual

df$times.rgc
df$times.ugc
df$ir100.sti.prep

plot(sim, y = "i.prev")
plot(sim, y = c("prev.rgcct", "prev.ugcct"),
     mean.col = 1:2, leg = TRUE)
plot(sim, y = c("ir100.gc", "ir100.ct"))


# Testing/Timing ------------------------------------------------------

control$bi.mods

undebug(sti_recov)
undebug(sti_tx)
undebug(prevalence_msm)
undebug(sti_trans)

dat <- initialize_msm(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- aging_msm(dat, at)       ## <1 ms
  dat <- deaths_msm(dat, at)      ## 4 ms
  dat <- births_msm(dat, at)      ## 6 ms
  dat <- test_msm(dat, at)        ## 2 ms
  dat <- tx_msm(dat, at)          ## 3 ms
  dat <- prep_msm(dat, at)        ## 2 ms
  dat <- progress_msm(dat, at)    ## 2 ms
  dat <- vl_msm(dat, at)          ## 3 ms
  dat <- simnet_msm(dat, at)      ## 53 ms
  dat <- disclose_msm(dat, at)    ## 1 ms
  dat <- acts_msm(dat, at)        ## 1 ms
  dat <- condoms_msm(dat, at)     ## 2 ms
  dat <- riskhist_msm(dat, at)    ## 4 ms
  dat <- position_msm(dat, at)    ## 1 ms
  dat <- trans_msm(dat, at)       ## 1 ms
  dat <- sti_trans(dat, at)       ## 4 ms
  dat <- sti_recov(dat, at)       ## 3 ms
  dat <- sti_tx(dat, at)          ## 2 ms
  dat <- prevalence_msm(dat, at)  ## 1 ms
  cat(at, ".", sep = "")
}

library(microbenchmark)

res <- microbenchmark(simnet_msm(dat, at = 2), times = 100)
summary(res, unit = "ms")

