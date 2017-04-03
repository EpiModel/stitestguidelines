
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
                     syph.rhiv.rr = 1.00,
                     syph.uhiv.rr = 1.00,
                     
                     rgc.tprob = 0.41333000,
                     ugc.tprob = 0.31404720,
                     rct.tprob = 0.19075540,
                     uct.tprob = 0.16394697,
                     
                     hivdx.syph.sympt.tx.rr = 1.45,
                     
                     hiv.rgc.rr = x[6],
                     hiv.ugc.rr = x[7],
                     hiv.rct.rr = x[6],
                     hiv.uct.rr = x[7],
                     
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
  prev.primsecosyph <- mean(df$prev.primsecosyph)
  #prev.primsecosyph.hivpos <- mean(df$prev.primsecosyph.hivpos)
  #prev.primsecosyph.hivneg <- mean(df$prev.primsecosyph.hivneg)
  #prev.hiv.primsecosyphpos <- mean(df$prev.hiv.primsecosyphpos)
  gcslope <- mean(df$ir100.gc[52] - df$ir100.gc[42])
  ctslope <- mean(df$ir100.ct[52] - df$ir100.ct[42])
  hivslope <- mean(df$ir100[52] - df$ir100[42])
  hivprevslope <- mean(df$i.prev[52] - df$i.prev[42])
  syphprevslope <- mean(df$prev.primsecosyph[52] - df$prev.primsecosyph[42])

  out <- c(gc.incid, ct.incid, hiv.incid, hiv.prev, prev.primsecosyph,
           #prev.primsecosyph.hivpos, prev.primsecosyph.hivneg, prev.hiv.primsecosyphpos,
           gcslope, ctslope, hivslope, hivprevslope, syphprevslope)

  return(out)
}


priors <- list(c("unif", 0.040, 0.060),
               c("unif", 0.030, 0.050),
               c("unif", 2.80, 3.10),
               c("unif", 1.50, 2.50),
               c("unif", 2.20, 2.80),
               c("unif", 1.20, 1.80))

targets <- c(4.2, 6.6, 3.8, 0.26, 0.046, 0, 0, 0, 0, 0) #0.103, 0.026,0.498,


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
