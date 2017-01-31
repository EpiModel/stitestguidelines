
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
                     
                     ai.scale = x[2],
                     rsyph.tprob = x[3],
                     usyph.tprob = x[4],
                     hiv.syph.rr = x[5],
                     syph.hiv.rr = x[6],
                     rgc.tprob = x[7],
                     ugc.tprob = x[8],
                     rct.tprob = x[9],
                     uct.tprob = x[10]
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
  syph.incid <- mean(df$ir100.syph)
  hiv.prev <- mean(df$i.prev)
  prev.syph.hivpos <- mean(df$prev.syph.hivpos)
  prev.syph.hivneg <- mean(df$prev.syph.hivneg)
  prev.syph <- mean(df$prev.syph)

  out <- c(gc.incid, ct.incid, hiv.incid, syph.incid, hiv.prev, prev.syph.hivpos, prev.syph.hivneg, prev.syph)

  return(out)
}


priors <- list(c("unif", 1.115, 1.130),
               c("unif", 0.020, 0.030),
               c("unif", 0.010, 0.020),
               c("unif", 1.90, 2.30),
               c("unif", 1.90, 2.30),
               c("unif", 0.37, 0.41),
               c("unif", 0.23, 0.28),
               c("unif", 0.27, 0.33),
               c("unif", 0.22, 0.28))

targets <- c(4.2, 6.6, 3.8, 0.9, 0.26, 0.103, 0.026, 0.046)


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
