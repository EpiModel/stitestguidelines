
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
                     
                     #ai.scale = x[2],
                     rsyph.tprob = x[2],
                     usyph.tprob = x[3],
                     hiv.rsyph.rr = x[4],
                     hiv.usyph.rr = x[5],
                     syph.hiv.rr = x[6],
                     rgc.tprob = x[7],
                     ugc.tprob = x[8],
                     rct.tprob = x[9],
                     uct.tprob = x[10],
                     hiv.rct.rr = x[4],
                     hiv.uct.rr = x[5],
                     hiv.rgc.rr = x[4],
                     hiv.ugc.rr = x[5],
                     syph.prim.sympt.prob.tx = x[11],
                     syph.seco.sympt.prob.tx = x[12],
                     syph.earlat.prob.tx = x[13],
                     syph.latelat.prob.tx = x[14]
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
  #syph.incid <- mean(df$ir100.syph)
  hiv.prev <- mean(df$i.prev)
  prev.syph.hivpos <- mean(df$prev.syph.hivpos)
  prev.syph.hivneg <- mean(df$prev.syph.hivneg)
  prev.syph <- mean(df$prev.syph)
  prev.hiv.syphpos <- mean(df$prev.hiv.syphpos)
  prev.earlysyph <- mean(df$prev.earlysyph)
  prev.latesyph <- mean(df$prev.latesyph)
  # prev.stage.incubprim <- mean(df$prev.stage.incubprim)
  # prev.stage.seco <- mean(df$prev.stage.seco)
  # prev.stage.earlat <- mean(df$prev.stage.earlat)
  # prev.stage.latelat <- mean(df$prev.stage.latelat)
  # prev.stage.latelatelat <- mean(df$prev.stage.latelatelat)
  # prev.stage.tert <- mean(df$prev.stage.tert)
  # prev.earlysyph <- mean(df$prev.earlysyph)
  # prev.latesyph <- mean(df$prev.latesyph)

  out <- c(gc.incid, ct.incid, hiv.incid, #syph.incid,
           hiv.prev, prev.syph.hivpos, prev.syph.hivneg, prev.syph,
           prev.hiv.syphpos, prev.earlysyph, prev.latesyph)
           #prev.stage.incubprim, prev.stage.seco, prev.stage.earlat, prev.stage.latelat, prev.stage.latelatelat,
           #prev.stage.tert)

  return(out)
}


priors <- list(#c("unif", 1.118, 1.130),
               c("unif", 0.0465, 0.050),
               c("unif", 0.0395, 0.042),
               c("unif", 2.30, 2.80),
               c("unif", 1.505, 1.78),
               c("unif", 1.70, 1.91),
               c("unif", 0.35, 0.415),
               c("unif", 0.28, 0.36),
               c("unif", 0.165, 0.220),
               c("unif", 0.167, 0.197),
               c("unif", 0.355, 0.415),
               c("unif", 0.625, 0.80),
               c("unif", 0.10, 0.25),
               c("unif", 0.08, 0.185))

targets <- c(4.2, 6.6, 3.8, #0.9,
             0.26, 0.103, 0.026, 0.046, 0.498, 0.554, 0.446)# #0.1385, 0.1385, 0.277, 0.20, 0.20, 0.046)


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
