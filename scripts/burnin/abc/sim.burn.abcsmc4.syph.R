
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
                     syph.tprob = x[3],
                     hiv.syph.rr = x[4],
                     syph.hiv.rr = x[5],
                     rgc.tprob = x[6],
                     ugc.tprob = x[7],
                     rct.tprob = x[8],
                     uct.tprob = x[9]
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

  out <- c(gc.incid, ct.incid, hiv.incid, syph.incid, hiv.prev)

  return(out)
}


priors <- list(c("unif", 1.118, 1.120),
               c("unif", 0.019, 0.020),
               c("unif", 1.90, 2.20),
               c("unif", 2.10, 2.30),
               c("unif", 0.37, 0.41),
               c("unif", 0.23, 0.28),
               c("unif", 0.27, 0.33),
               c("unif", 0.22, 0.28))

targets <- c(4.2, 6.6, 3.8, 0.9, 0.26)


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
