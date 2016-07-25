
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
    rgc.tprob <- runif(1, 0.35, 0.60)
    ugc.tprob <- runif(1, 0.20, 0.40)
    rct.tprob <- runif(1, 0.35, 0.60)
    uct.tprob <- runif(1, 0.20, 0.40)

    rgc.dur.asympt <- runif(1, 26, 52)
    ugc.dur.asympt <- runif(1, 26, 52)
    rct.dur.asympt <- runif(1, 39, 65)
    uct.dur.asympt <- runif(1, 39, 65)

    rgc.sympt.prob <- runif(1, 0.05, 0.20)
    ugc.sympt.prob <- runif(1, 0.60, 0.95)
    rct.sympt.prob <- runif(1, 0.05, 0.20)
    uct.sympt.prob <- runif(1, 0.60, 0.95)

    hiv.rect.rr <- runif(1, 2, 3)
    hiv.ureth.rr <- runif(1, 1, 2)
  }
  if (batch > 1) {
    load(paste0("simChosen.b", batch-1, ".rda"))
    mn <- apply(simChosen, 2, mean)
    sds <- apply(simChosen, 2, sd)
    lo <- mn - sds
    hi <- mn + sds

    rgc.tprob <- runif(1, lo[["rgc.tprob"]], hi[["rgc.tprob"]])
    ugc.tprob <- runif(1, lo[["ugc.tprob"]], hi[["ugc.tprob"]])
    rct.tprob <- runif(1, lo[["rct.tprob"]], hi[["rct.tprob"]])
    uct.tprob <- runif(1, lo[["uct.tprob"]], hi[["uct.tprob"]])

    rgc.dur.asympt <- runif(1, lo[["rgc.dur.asympt"]], hi[["rgc.dur.asympt"]])
    ugc.dur.asympt <- runif(1, lo[["ugc.dur.asympt"]], hi[["ugc.dur.asympt"]])
    rct.dur.asympt <- runif(1, lo[["rct.dur.asympt"]], hi[["rct.dur.asympt"]])
    uct.dur.asympt <- runif(1, lo[["uct.dur.asympt"]], hi[["uct.dur.asympt"]])

    rgc.sympt.prob <- runif(1, lo[["rgc.sympt.prob"]], hi[["rgc.sympt.prob"]])
    ugc.sympt.prob <- runif(1, lo[["ugc.sympt.prob"]], hi[["ugc.sympt.prob"]])
    rct.sympt.prob <- runif(1, lo[["rct.sympt.prob"]], hi[["rgc.sympt.prob"]])
    uct.sympt.prob <- runif(1, lo[["uct.sympt.prob"]], hi[["uct.sympt.prob"]])

    hiv.rect.rr <- runif(1, lo[["hiv.rect.rr"]], hi[["hiv.rect.rr"]])
    hiv.ureth.rr <- runif(1, lo[["hiv.ureth.rr"]], hi[["hiv.ureth.rr"]])
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
                     gc.dur.tx = 2,
                     gc.dur.ntx = NULL,

                     rct.dur.asympt = rct.dur.asympt,
                     uct.dur.asympt = uct.dur.asympt,
                     ct.dur.tx = 2,
                     ct.dur.ntx = NULL,

                     gc.prob.cease = 0,
                     ct.prob.cease = 0,

                     gc.sympt.prob.tx = 0.90,
                     ct.sympt.prob.tx = 0.85,
                     gc.asympt.prob.tx = 0,
                     ct.asympt.prob.tx = 0,

                     prep.sti.screen.int = 182,
                     prep.sti.prob.tx = 1,

                     sti.cond.rr = 0.3,

                     hiv.rgc.rr = hiv.rect.rr,
                     hiv.ugc.rr = hiv.ureth.rr,
                     hiv.rct.rr = hiv.rect.rr,
                     hiv.uct.rr = hiv.ureth.rr,
                     hiv.dual.rr = 0)

  init <- init_msm(nwstats = st,
                   prev.B = 0.253,
                   prev.W = 0.253,
                   prev.ugc = 0.10,
                   prev.rgc = 0.10,
                   prev.uct = 0.10,
                   prev.rct = 0.10)

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
  rect.prev <- mean(df$prev.rgcct)
  ureth.prev <- mean(df$prev.ugcct)
  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.prev <- mean(df$i.prev)

  out <- data.frame(rgc.tprob = rgc.tprob,
                    ugc.tprob = ugc.tprob,
                    rct.tprob = rct.tprob,
                    uct.tprob = uct.tprob,
                    rgc.dur.asympt = rgc.dur.asympt,
                    ugc.dur.asympt = ugc.dur.asympt,
                    rct.dur.asympt = rct.dur.asympt,
                    uct.dur.asympt = uct.dur.asympt,
                    rgc.sympt.prob = rgc.sympt.prob,
                    ugc.sympt.prob = ugc.sympt.prob,
                    rct.sympt.prob = rct.sympt.prob,
                    uct.sympt.prob = uct.sympt.prob,
                    hiv.rect.rr = hiv.rect.rr,
                    hiv.ureth.rr = hiv.ureth.rr,
                    rect.prev = rect.prev,
                    ureth.prev = ureth.prev,
                    gc.incid = gc.incid,
                    ct.incid = ct.incid,
                    hiv.prev = hiv.prev)

  return(out)
}

rejection <- function(sim, targets, threshold) {
  edist <- sapply(1:nrow(sim), function(x) sqrt(sum((targets - sim[x, 15:19])^2)))
  edist.quant <- quantile(edist, threshold)
  accepted <- which(edist <= edist.quant)
  post <- sim[accepted, ]
  cat("\n Accepted n:", nrow(post))
  return(post)
}


# Parameters ----------------------------------------------------------

target.n.chosen <- NULL

# rect.prev, ureth.prev, gc.incid, ct.incid, hiv.prev
targets <- c(0.17, 0.07, 43, 48, 0.26)

threshold <- 0.01
sims.per.batch <- 100

batch <- 1
out.fn.all <- paste0("data/simDataAll.b", batch, ".rda")
out.fn.chosen <- paste0("data/simDataChosen.b", batch, ".rda")


# Algorithm -----------------------------------------------------------

# if target.n.chosen is NULL, then do the rejection manually
if (is.null(target.n.chosen)) {

  # Run batches of sims
  cl <- parallel::makeCluster(parallel::detectCores())
  registerDoParallel(cl)
  nsims <- sims.per.batch
  sout <- foreach(s = 1:nsims) %dopar% {
    f(batch = batch)
  }
  parallel::stopCluster(cl)
  sim <- as.data.frame(do.call("rbind", sout))

  ind.fn <- tempfile(pattern = paste0("simDataAll.b", batch, "."),
                     tmpdir = "data/ind/", fileext = ".rda")

  # Save all sims so far
  if (!file.exists(out.fn.all)) {
    save(sim, file = out.fn.all)
    save(sim, file = ind.fn)
  } else {
    save(sim, file = ind.fn)
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
    nsims <- sims.per.batch
    sout <- foreach(s = 1:nsims) %dopar% {
      f(batch = batch)
    }
    stopCluster(cl)
    sim <- as.data.frame(do.call("rbind", sout))

    ind.fn <- tempfile(pattern = paste0("simDataAll.b", batch, "."),
                       tmpdir = "data/ind/", fileext = ".rda")

    # Rejection algorithm
    simChosen <- rejection(sim, targets = targets, threshold = threshold)

    # Save all sims so far
    if (!file.exists(out.fn.all)) {
      save(sim, file = out.fn.all)
      save(sim, file = ind.fn)
    } else {
      save(sim, file = ind.fn)
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
