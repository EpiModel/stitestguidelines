
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
                     rgc.tprob = x[8],
                     ugc.tprob = x[9],
                     rct.tprob = x[10],
                     uct.tprob = x[11],
                     hivdx.syph.sympt.tx.rr = x[12]
                     # hiv.rct.rr = x[4],
                     # hiv.uct.rr = x[5],
                     # hiv.rgc.rr = x[4],
                     # hiv.ugc.rr = x[5],
                     # syph.prim.sympt.prob.tx = x[12],
                     # syph.seco.sympt.prob.tx = x[13],
                     # syph.earlat.prob.tx = x[14],
                     # syph.latelat.prob.tx = x[15]
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
  prev.primsecosyph.hivpos <- mean(df$prev.primsecosyph.hivpos)
  prev.primsecosyph.hivneg <- mean(df$prev.primsecosyph.hivneg)
  prev.primsecosyph <- mean(df$prev.primsecosyph)
  prev.hiv.primsecosyphpos <- mean(df$prev.hiv.primsecosyphpos)
  #prev.earlysyph <- mean(df$prev.earlysyph)
  #prev.latesyph <- mean(df$prev.latesyph)

  out <- c(gc.incid, ct.incid, hiv.incid, #syph.incid,
           hiv.prev, prev.primsecosyph.hivpos, prev.primsecosyph.hivneg, prev.primsecosyph,
           prev.hiv.primsecosyphpos)#, prev.earlysyph, prev.latesyph)

  return(out)
}


priors <- list(#c("unif", 1.118, 1.130),
               c("unif", 0.050, 0.080),
               c("unif", 0.040, 0.070),
               c("unif", 2.50, 3.50),
               c("unif", 1.50, 2.50),
               c("unif", 4.00, 10.00),
               c("unif", 2.00, 8.00),
               c("unif", 0.35, 0.50),
               c("unif", 0.25, 0.40),
               c("unif", 0.15, 0.22),
               c("unif", 0.12, 0.20),
               c("unif", 1.70, 3.00))
               # c("unif", 0.30, 0.60),
               # c("unif", 0.60, 0.80),
               # c("unif", 0.10, 0.25),
               # c("unif", 0.05, 0.20))

targets <- c(4.2, 6.6, 3.8, #0.9,
             0.26, 0.103, 0.026, 0.046, 0.498)#, 0.554, 0.446)# #0.1385, 0.1385, 0.277, 0.20, 0.20, 0.046)


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
