
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))

f <- function() {

  sourceDir("source/", verbose = FALSE)

  rgc.tprob = runif(1, 0.2, 0.5)
  ugc.tprob = runif(1, 0.2, 0.5)
  rct.tprob = runif(1, 0.2, 0.5)
  uct.tprob = runif(1, 0.2, 0.5)

  rgc.dur.asympt = runif(1, 4*4.33, 12*4.33)
  ugc.dur.asympt = runif(1, 4*4.33, 12*4.33)
  gc.dur.ntx = runif(1, 4*4.33, 12*4.33)

  rct.dur.asympt = runif(1, 4*4.33, 12*4.33)
  uct.dur.asympt = runif(1, 4*4.33, 12*4.33)
  ct.dur.ntx = runif(1, 4*4.33, 12*4.33)

  rgc.sympt.prob = runif(1, 0.05, 0.25)
  ugc.sympt.prob = runif(1, 0.5, 0.9)
  rct.sympt.prob = runif(1, 0.05, 0.25)
  uct.sympt.prob = runif(1, 0.5, 0.9)

  load("est/nwstats.rda")
  param <- param_msm(nwstats = st,
                     ai.scale = 1,

                     prep.coverage = 0,

                     rcomp.prob = 0,
                     rcomp.adh.groups = 0:3,
                     rcomp.main.only = FALSE,
                     rcomp.discl.only = FALSE,

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
                     gc.dur.tx = 14/7,
                     gc.dur.ntx = gc.dur.ntx,

                     rct.dur.asympt = rct.dur.asympt,
                     uct.dur.asympt = uct.dur.asympt,
                     ct.dur.tx = 14/7,
                     ct.dur.ntx = ct.dur.ntx,

                     gc.prob.cease = 0,
                     ct.prob.cease = 0,

                     gc.sympt.prob.tx = 0.90,
                     ct.sympt.prob.tx = 0.85,
                     gc.asympt.prob.tx = 0,
                     ct.asympt.prob.tx = 0,

                     prep.sti.screen.int = 182,
                     prep.sti.prob.tx = 1,

                     sti.cond.rr = 0.3,

                     hiv.rgc.rr = 2.5,
                     hiv.ugc.rr = 1.25,
                     hiv.rct.rr = 2.5,
                     hiv.uct.rr = 1.25,
                     hiv.dual.rr = 0)

  init <- init_msm(nwstats = st,
                   prev.B = 0.253,
                   prev.W = 0.253,
                   prev.ugc = 0.1,
                   prev.rgc = 0.1,
                   prev.uct = 0.1,
                   prev.rct = 0.1)

  control <- control_msm(simno = 1,
                         nsteps = 1300,
                         nsims = 1,
                         ncores = 1,
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


  load("est/fit.rda")
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 500)
  rgc.prev <- mean(df$prev.rgc)
  ugc.prev <- mean(df$prev.ugc)
  rct.prev <- mean(df$prev.rct)
  uct.prev <- mean(df$prev.uct)

  out <- data.frame(rgc.tprob = rgc.tprob, ugc.tprob = ugc.tprob,
                    rct.tprob = rct.tprob, uct.tprob = uct.tprob,
                    rgc.dur.asympt = rgc.dur.asympt, ugc.dur.asympt = ugc.dur.asympt, gc.dur.ntx = gc.dur.ntx,
                    rct.dur.asympt = rct.dur.asympt, uct.dur.asympt = uct.dur.asympt, ct.dur.ntx = ct.dur.ntx,
                    rgc.sympt.prob = rgc.sympt.prob, ugc.sympt.prob = ugc.sympt.prob,
                    rct.sympt.prob = rct.sympt.prob, uct.sympt.prob = uct.sympt.prob,
                    rgc.prev = rgc.prev, ugc.prev = ugc.prev,
                    rct.prev = rct.prev, uct.prev = uct.prev)

  return(out)
}

# rejection <- function(sim, target.stat = 0.26, threshold = 0.0005) {
#   edist <- sapply(1:nrow(sim), function(x) sqrt(sum((target.stat - sim$stat.mean[x])^2)))
#   accepted <- which(edist <= threshold)
#   post <- sim[accepted, ]
#   return(post)
# }

out.fn.all <- "data/simDataAll.rda"
# out.fn.chosen <- "data/simDataChosen.rda"

# if (file.exists(out.fn.chosen)) {
#   while(inherits(try(load(out.fn.chosen), silent = TRUE), "try-error")) Sys.sleep(1)
#   n.chosen <- nrow(simChosen)
# } else {s
#   n.chosen <- 0
# }

# target.n.chosen <- 100

# while (n.chosen < target.n.chosen) {

  registerDoParallel(parallel::detectCores())
  nsims <- 50
  sout <- foreach(s = 1:nsims) %dopar% {
    f()
  }

  sim <- as.data.frame(do.call("rbind", sout))
  # simChosen <- rejection(sim)

  # Save all sims
  if (!file.exists(out.fn.all)) {
    save(sim, file = out.fn.all)
  } else {
    simNew <- sim
    while(inherits(try(load(out.fn.all), silent = TRUE), "try-error")) Sys.sleep(1)
    sim <- rbind(sim, simNew)
    save(sim, file = out.fn.all)
  }

  # tot.sim <- nrow(sim)

  # Save accepted sims
  # if (!file.exists(out.fn.chosen) & nrow(simChosen) > 0) {
  #   save(simChosen, file = out.fn.chosen)
  # }
  # if (file.exists(out.fn.chosen)) {
  #   if (nrow(simChosen) == 0) {
  #     load(out.fn.chosen)
  #   } else {
  #     simChosenNew <- simChosen
  #     while(inherits(try(load(out.fn.chosen), silent = TRUE), "try-error")) Sys.sleep(1)
  #     simChosen <- rbind(simChosen, simChosenNew)
  #     save(simChosen, file = out.fn.chosen)
  #     if (nrow(simChosen) > target.n.chosen) {
  #       samp <- sample(1:nrow(simChosen), size = target.n.chosen)
  #       simChosen <- simChosen[samp, , drop = FALSE]
  #     }
  #   }
  # }

  # n.chosen <- nrow(simChosen)
  # p.chosen <- round(n.chosen / tot.sim, 3)
  # cat("\n tot.sim=", tot.sim, " n.chosen=", n.chosen, " p.chosen=", p.chosen, sep = "")

# }

