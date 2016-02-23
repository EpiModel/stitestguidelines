
## Test Script for stiPrEP Project

library(EpiModelHIVmsm)
sourceDir("source/", TRUE)

load("est/nwstats.10k.rda")

param <- param.msm(nwstats = st,
                   testing.pattern = "interval",
                   ai.scale = 1.323,
                   riskh.start = 1,
                   prep.start = 30,
                   prep.elig.model = "cdc3",
                   prep.class.prob = reallocate_pcp(reall = 0),
                   prep.class.hr = c(1, 0.69, 0.19, 0.05),
                   prep.coverage = 0.5,
                   prep.cov.method = "curr",
                   prep.cov.rate = 1,
                   prep.tst.int = 90,
                   prep.risk.int = 182,
                   rcomp.prob = 0.5,
                   rcomp.hadhr.only = FALSE,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE)
init <- init.msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253)
control <- control.msm(simno = 1,
                       nsteps = 100,
                       nsims = 1,
                       ncores = 1,
                       save.int = 5000,
                       verbose.int = 1,
                       save.other = c("attr", "temp", "riskh", "el", "p"),
                       acts.FUN = acts.sti,
                       condoms.FUN = condoms.sti,
                       initialize.FUN = initialize.sti,
                       prep.FUN = prep.sti,
                       prev.FUN = prevalence.sti,
                       riskhist.FUN = riskhist.msm,
                       trans.FUN = trans.msm,
                       test.FUN = test.sti)

load("est/fit.10k.rda")
sim <- netsim(est, param, init, control)

