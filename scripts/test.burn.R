
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
# devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

# Main Test Script ----------------------------------------------------

load("est/nwstats.rda")

load("est/abc.avg.parms.1pct.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,
                   prep.start = 2601)
init <- init_msm(nwstats = st)

control <- control_msm(simno = 1,
                       nsteps = 2600,
                       nsims = 8,
                       ncores = 8)

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


plot(sim, y = "i.prev")
plot(sim, y = c("prev.rgcct", "prev.ugcct"),
     mean.col = 1:2, leg = TRUE)
plot(sim, y = c("ir100.gc", "ir100.ct"))


# Testing/Timing ------------------------------------------------------

control$bi.mods

debug(sti_recov)
debug(sti_tx)

dat <- initialize_msm(est, param, init, control, s = 1)

for (at in 2:30) {
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

