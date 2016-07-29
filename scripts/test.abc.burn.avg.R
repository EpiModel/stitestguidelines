suppressMessages(library("EpiModelHIV"))
sourceDir("source/", verbose = FALSE)

param <- param_msm(nwstats = st,
                   ai.scale = 1,

                   prep.coverage = 0,
                   prep.start = 1e8,

                   rcomp.prob = 0,
                   rcomp.adh.groups = 0:3,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE,

                   rgc.tprob = 0.43,
                   ugc.tprob = 0.20,
                   rct.tprob = 0.36,
                   uct.tprob = 0.17,

                   rgc.sympt.prob = 0.09137000 ,
                   ugc.sympt.prob = 0.71578944,
                   rct.sympt.prob = 0.06695481,
                   uct.sympt.prob = 0.88658264,

                   rgc.dur.asympt = 35.12593426 ,
                   ugc.dur.asympt = 35.67146136 ,
                   gc.dur.tx = 2,
                   gc.dur.ntx = NULL,

                   rct.dur.asympt = 52.33361331,
                   uct.dur.asympt = 51.65952974,
                   ct.dur.tx = 2,
                   ct.dur.ntx = NULL,

                   gc.prob.cease = 0.11,
                   ct.prob.cease = 0.09,

                   gc.sympt.prob.tx = 0.90,
                   ct.sympt.prob.tx = 0.90,
                   gc.asympt.prob.tx = 0,
                   ct.asympt.prob.tx = 0,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,
                   prep.continue.stand.tx = TRUE,

                   sti.cond.rr = 0.3,

                   hiv.rgc.rr = 2.8,
                   hiv.ugc.rr = 1.8,
                   hiv.rct.rr = 2.8,
                   hiv.uct.rr = 1.8,
                   hiv.dual.rr = 0.2)

init <- init_msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.02,
                 prev.rgc = 0.05,
                 prev.uct = 0.02,
                 prev.rct = 0.05)

control <- control_msm(simno = 1,
                       nsteps = 1200,
                       nsims = 16,
                       ncores = 16,
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
                       verbose = FALSE,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN",
                                        "test.FUN", "tx.FUN", "prep.FUN",
                                        "progress.FUN", "vl.FUN",
                                        "resim_nets.FUN", "disclose.FUN",
                                        "acts.FUN", "condoms.FUN", "riskhist.FUN",
                                        "position.FUN", "trans.FUN", "stitrans.FUN",
                                        "stirecov.FUN", "stitx.FUN", "prev.FUN"))


data(est)
sim <- netsim(est, param, init, control)


df <- tail(as.data.frame(sim), 200)
rect.prev <- mean(df$prev.rgcct)
ureth.prev <- mean(df$prev.ugcct)
gc.incid <- mean(df$ir100.gc)
ct.incid <- mean(df$ir100.ct)
hiv.incid <- mean(df$ir100)
hiv.prev <- mean(df$i.prev)

out <- c(rect.prev, ureth.prev, gc.incid, ct.incid, hiv.incid, hiv.prev)
names(out) <- c("rect.prev", "ureth.prev", "gc.incid", "ct.incid", "hiv.incid", "hiv.prev")
targets <- c(0.135, 0.046, 23.2, 26.8, 3.8, 0.26)
data.frame(out, targets)

plot(sim, y = c("prev.rgcct", "prev.ugcct"),
     mean.col = 1:2, qnts = FALSE, leg = TRUE)
plot(sim, y = c("ir100.gc", "ir100.ct"),
     mean.col = 1:2, qnts = FALSE, leg = TRUE)
plot(sim, y = "i.prev")
