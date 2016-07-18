
# Packages ------------------------------------------------------------

library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))


# Functions -----------------------------------------------------------


f <- function(batch) {

  suppressMessages(library("EpiModelHIV"))
  sourceDir("source/", verbose = FALSE)

  if (batch == 1) {
    rgc.tprob <- runif(1, 0.2, 0.6)
    ugc.tprob <- runif(1, 0.2, 0.5)
    rct.tprob <-runif(1, 0.2, 0.6)
    uct.tprob <- runif(1, 0.2, 0.5)

    rgc.dur.asympt <- runif(1, 17, 52)
    ugc.dur.asympt <- runif(1, 17, 52)
    gc.dur.ntx <- mean(rgc.dur.asympt, ugc.dur.asympt)

    rct.dur.asympt <- runif(1, 17, 52)
    uct.dur.asympt <- runif(1, 17, 52)
    ct.dur.ntx <- mean(rgc.dur.asympt, ugc.dur.asympt)

    rgc.sympt.prob <- runif(1, 0.05, 0.25)
    ugc.sympt.prob <- runif(1, 0.50, 0.90)
    rct.sympt.prob <- runif(1, 0.05, 0.25)
    uct.sympt.prob <- runif(1, 0.50, 0.90)
  }
  if (batch > 1) {
    load(paste0("simChosen.b", batch-1, ".rda"))
    mn <- apply(simChosen, 2, mean)
    sds <- apply(simChosen, 2, sd)
    lo <- mn - sds
    hi <- mn + sds

    rgc.tprob <- runif(1, lo[["rgc.tprob"]], hi[["rgc.tprob"]])
    ugc.tprob <- runif(1, lo[["ugc.tprob"]], hi[["ugc.tprob"]])
    rct.tprob <-runif(1, lo[["rct.tprob"]], hi[["rct.tprob"]])
    uct.tprob <- runif(1, lo[["uct.tprob"]], hi[["uct.tprob"]])

    rgc.dur.asympt <- runif(1, lo[["rgc.dur.asympt"]], hi[["rgc.dur.asympt"]])
    ugc.dur.asympt <- runif(1, lo[["ugc.dur.asympt"]], hi[["ugc.dur.asympt"]])
    gc.dur.ntx <- mean(rgc.dur.asympt, ugc.dur.asympt)

    rct.dur.asympt <- runif(1, lo[["rct.dur.asympt"]], hi[["rct.dur.asympt"]])
    uct.dur.asympt <- runif(1, lo[["uct.dur.asympt"]], hi[["uct.dur.asympt"]])
    ct.dur.ntx <- mean(rgc.dur.asympt, ugc.dur.asympt)

    rgc.sympt.prob <- runif(1, lo[["rgc.sympt.prob"]], hi[["rgc.sympt.prob"]])
    ugc.sympt.prob <- runif(1, lo[["ugc.sympt.prob"]], hi[["ugc.sympt.prob"]])
    rct.sympt.prob <- runif(1, lo[["rct.sympt.prob"]], hi[["rgc.sympt.prob"]])
    uct.sympt.prob <- runif(1, lo[["uct.sympt.prob"]], hi[["uct.sympt.prob"]])
  }

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
                   prev.ugc = 0.111,
                   prev.rgc = 0.102,
                   prev.uct = 0.084,
                   prev.rct = 0.141)

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
                    rgc.dur.asympt = rgc.dur.asympt, ugc.dur.asympt = ugc.dur.asympt,
                    rct.dur.asympt = rct.dur.asympt, uct.dur.asympt = uct.dur.asympt,
                    rgc.sympt.prob = rgc.sympt.prob, ugc.sympt.prob = ugc.sympt.prob,
                    rct.sympt.prob = rct.sympt.prob, uct.sympt.prob = uct.sympt.prob,
                    rgc.prev = rgc.prev, ugc.prev = ugc.prev,
                    rct.prev = rct.prev, uct.prev = uct.prev)

  return(out)
}

rejection <- function(sim, targets, threshold) {
  diff.rgc <- abs(sim$rgc.prev - targets[1])
  diff.ugc <- abs(sim$ugc.prev - targets[2])
  diff.rct <- abs(sim$rct.prev - targets[3])
  diff.uct <- abs(sim$uct.prev - targets[4])

  choice <- which(diff.rgc <= threshold & diff.ugc <= threshold &
                    diff.rct <= threshold & diff.uct <= threshold)
  simChosen <- sim[choice, , drop = FALSE]
  return(simChosen)
}


# Parameters ----------------------------------------------------------

target.n.chosen <- 100
targets <- c(0.102, 0.111, 0.141, 0.084)
threshold <- 0.002
sims.per.batch <- 16
batch <- 4
out.fn.all <- paste0("data/simDataAll.b", batch, ".rda")
out.fn.chosen <- paste0("data/simDataChosen.b", batch, ".rda")



# Algorithm -----------------------------------------------------------

# if target.n.chosen is NULL, then do the rejection manually
if (is.null(target.n.chosen)) {

  # Run batches of sims
  cl <- makeCluster(parallel::detectCores())
  registerDoParallel(cl)
  nsims <- sims.per.batch
  sout <- foreach(s = 1:nsims) %dopar% {
    f(batch = batch)
  }
  stopCluster(cl)
  sim <- as.data.frame(do.call("rbind", sout))

  # Save all sims so far
  if (!file.exists(out.fn.all)) {
    save(sim, file = out.fn.all)
  } else {
    simNew <- sim
    while(inherits(try(load(out.fn.all), silent = TRUE), "try-error")) {
      Sys.sleep(1)
    }
    sim <- rbind(sim, simNew)
    save(sim, file = out.fn.all)
  }

} else {

  # Load current simChosen file to get number already chosen
  if (file.exists(out.fn.chosen)) {
    while(inherits(try(load(out.fn.chosen), silent = TRUE), "try-error")) {
      Sys.sleep(1)
    }
    n.chosen <- nrow(simChosen)
  } else {
    n.chosen <- 0
  }

  # ABC-R Loop
  while (n.chosen < target.n.chosen) {

    # Run batches of sims
    cl <- makeCluster(parallel::detectCores())
    registerDoParallel(cl)
    nsims <- parallel::detectCores()
    sout <- foreach(s = 1:nsims) %dopar% {
      f(batch = batch)
    }
    stopCluster(cl)
    sim <- as.data.frame(do.call("rbind", sout))

    # Rejection algorithm
    simChosen <- rejection(sim, targets = targets, threshold = threshold)

    # Save all sims so far
    if (!file.exists(out.fn.all)) {
      save(sim, file = out.fn.all)
    } else {
      simNew <- sim
      while(inherits(try(load(out.fn.all), silent = TRUE), "try-error")) {
        Sys.sleep(1)
      }
      sim <- rbind(sim, simNew)
      save(sim, file = out.fn.all)
    }

    # Save accepted sims
    if (!file.exists(out.fn.chosen)) {
      save(simChosen, file = out.fn.chosen)
    }
    if (file.exists(out.fn.chosen)) {
      if (nrow(simChosen) == 0) {
        load(out.fn.chosen)
      } else {
        simChosenNew <- simChosen
        while(inherits(try(load(out.fn.chosen), silent = TRUE), "try-error")) {
          Sys.sleep(1)
        }
        simChosen <- rbind(simChosen, simChosenNew)
        save(simChosen, file = out.fn.chosen)
      }
    }

    # Update n chosen within loop
    n.chosen <- nrow(simChosen)
  }

}
