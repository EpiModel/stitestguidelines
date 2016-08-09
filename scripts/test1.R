
rm(list=ls())
suppressMessages(library("EpiModelHIV"))

# devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

# Main Test Script ----------------------------------------------------

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
                   hiv.uct.rr = hiv.ureth.rr)

init <- init_msm(nwstats = st)

control <- control_msm(simno = 1,
                       nsteps = 250,
                       nsims = 1,
                       ncores = 1)

load("est/fit.rda")
sim <- netsim(est, param, init, control)

df <- as.data.frame(sim)
names(df)
df$incid.gc
df$ir100.gc
df$prev.gc.dual
df$prev.ct.dual

df$times.rgc
df$times.ugc

plot(sim, y = c("prev.rgc", "prev.ugc", "prev.rct", "prev.uct"),
     mean.col = 1:4, leg = TRUE)



# Testing/Timing ------------------------------------------------------

dat <- initialize_sti(est, param, init, control, s = 1)

for (at in 2:30) {
  dat <- aging_msm(dat, at)       ## <1 ms
  dat <- deaths_msm(dat, at)      ## 4 ms
  dat <- births_msm(dat, at)      ## 6 ms
  dat <- test_msm(dat, at)        ## 2 ms
  dat <- tx_msm(dat, at)          ## 3 ms
  dat <- prep_sti(dat, at)        ## 2 ms
  dat <- progress_msm(dat, at)    ## 2 ms
  dat <- vl_msm(dat, at)          ## 3 ms
  dat <- simnet_msm(dat, at)      ## 53 ms
  dat <- disclose_msm(dat, at)    ## 1 ms
  dat <- acts_sti(dat, at)        ## 1 ms
  dat <- condoms_sti(dat, at)     ## 2 ms
  dat <- riskhist_sti(dat, at)    ## 4 ms
  dat <- position_sti(dat, at)    ## 1 ms
  dat <- trans_sti(dat, at)       ## 1 ms
  dat <- sti_trans(dat, at)       ## 4 ms
  dat <- sti_recov(dat, at)       ## 3 ms
  dat <- sti_tx(dat, at)          ## 2 ms
  dat <- prevalence_msm(dat, at)  ## 1 ms
  cat(at, ".", sep = "")
}

library(microbenchmark)

res <- microbenchmark(simnet_msm(dat, at = 2), times = 100)
summary(res, unit = "ms")

