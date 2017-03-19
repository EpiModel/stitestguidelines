
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))
suppressMessages(library("EasyABC"))

f <- function(x) {

  set.seed(x[1])

  suppressMessages(library("EpiModelHIV"))

  data(st)
  param <- param_msm(nwstats = st,
                     
                     ai.scale = 1.11,
                     
                     rsyph.tprob = x[2],
                     usyph.tprob = x[3],
                     
                     hiv.rsyph.rr = x[4],
                     hiv.usyph.rr = x[5],
                     syph.rhiv.rr = x[6],
                     syph.uhiv.rr = x[7],
                     
                     rgc.tprob = 0.41333000,
                     ugc.tprob = 0.30904720,
                     rct.tprob = 0.19075540,
                     uct.tprob = 0.16394697,
                     
                     hivdx.syph.sympt.tx.rr = 1.45,
                     
                     hiv.rgc.rr = 2.35,
                     hiv.ugc.rr = 1.35,
                     hiv.rct.rr = 2.35,
                     hiv.uct.rr = 1.35,
                     
                     prep.coverage = 0,
                     stianntest.coverage = 0,
                     stihighrisktest.coverage = 0,
                     ept.coverage = 0,
                     
                     prep.start = 5000,
                     stitest.start = 5000,
                     ept.start = 5000
  )

  init <- init_msm(nwstats = st)

  control <- control_msm(simno = 1,
                         nsteps = 2600,
                         nsims = 1, ncores = 1,
                         verbose = FALSE)

  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 52)

  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.incid <- mean(df$ir100)
  hiv.prev <- mean(df$i.prev)
  prev.primsecosyph.hivpos <- mean(df$prev.primsecosyph.hivpos)
  prev.primsecosyph.hivneg <- mean(df$prev.primsecosyph.hivneg)
  prev.primsecosyph <- mean(df$prev.primsecosyph)
  prev.hiv.primsecosyphpos <- mean(df$prev.hiv.primsecosyphpos)
  gcslope <- mean(df$ir100.gc[52] - df$ir100.gc[32])
  ctslope <- mean(df$ir100.ct[52] - df$ir100.ct[32])
  hivslope <- mean(df$ir100[52] - df$ir100[32])
  hivprevslope <- mean(df$i.prev[52] - df$i.prev[32])
  syphprevslope <- mean(df$prev.primsecosyph[32] - df$prev.primsecosyph[32])

  out <- c(gc.incid, ct.incid, hiv.incid, hiv.prev, 
           prev.primsecosyph.hivpos, prev.primsecosyph.hivneg, prev.primsecosyph,
           prev.hiv.primsecosyphpos, gcslope, ctslope, hivslope, hivprevslope, syphprevslope)

  return(out)
}


priors <- list(c("unif", 0.040, 0.050),
               c("unif", 0.030, 0.040),
               c("unif", 2.50, 3.25),
               c("unif", 1.50, 2.00),
               c("unif", 6.00, 7.00),
               c("unif", 4.00, 5.50))

targets <- c(4.2, 6.6, 3.8, 0.26, 0.103, 0.026, 0.046, 0.498, 0, 0, 0, 0, 0)


( nsim <- as.numeric(Sys.getenv("NSIM")) )
( pacc <- as.numeric(Sys.getenv("PACC")) )

a <- ABC_sequential(method = "Lenormand",
                    model = f,
                    prior = priors,
                    nb_simul = nsim,
                    summary_stat_target = targets,
                    p_acc_min = pacc,
                    progress_bar = TRUE,
                    n_cluster = 16,
                    use_seed = TRUE,
                    verbose = FALSE)

fn <- paste0("data/smc4.", pacc*100, "pct.", nsim, "sim.rda")
save(a, file = fn)
