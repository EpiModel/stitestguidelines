
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

                     ai.scale = 1.04,

                     rsyph.tprob = x[2],
                     usyph.tprob = x[3],

                     hiv.rsyph.rr = x[4],
                     hiv.usyph.rr = x[5],
                     syph.rhiv.rr = 1.00,
                     syph.uhiv.rr = 1.00,

                     rgc.tprob = x[6],
                     ugc.tprob = x[7],
                     rct.tprob = x[8],
                     uct.tprob = x[9],

                     hiv.rgc.rr = x[4],
                     hiv.ugc.rr = x[5],
                     hiv.rct.rr = x[4],
                     hiv.uct.rr = x[5],

                     prep.coverage = 0,
                     stianntest.coverage = 0.1,
                     stihighrisktest.coverage = 0.1,
                     ept.coverage = 0,

                     prep.start = 5000,
                     stitest.start = 2601,
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
  hiv.prev <- mean(df$i.prev)
  syph.incid <- mean(df$ir100.syph)

  gcslope <- mean(df$ir100.gc[52] - df$ir100.gc[47])
  ctslope <- mean(df$ir100.ct[52] - df$ir100.ct[47])
  syphslope <- mean(df$ir100.syph[52] - df$ir100.syph[47])
  hivslope <- mean(df$ir100[52] - df$ir100[47])
  hivprevslope <- mean(df$i.prev[52] - df$i.prev[47])

  out <- c(gc.incid, ct.incid, hiv.prev, syph.incid,
           gcslope, ctslope, syphslope, hivslope,
           hivprevslope)


  return(out)
}


priors <- list(c("unif", 0.147, 0.149),
               c("unif", 0.127, 0.129),
               c("unif", 2.08, 2.12),
               c("unif", 1.28, 1.32),
               c("unif", 0.43, 0.45),
               c("unif", 0.33, 0.35),
               c("unif", 0.196, 0.20),
               c("unif", 0.176, 0.18))

targets <- c(3.5, 5.0, 0.15, 2.0, 0, 0, 0, 0, 0)

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
