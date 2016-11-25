
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

                     syph.tprob = x[2],
                     hiv.syph.rr = x[3],
                     syph.hiv.rr = x[4])

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
  syph.prev.hivpos <- mean(df$prev.syph.hivpos)
  syph.prev.hivneg <- mean(df$prev.syph.hivneg)
  syphratio <- (syph.prev.hivpos / syph.prev.hivneg)

  out <- c(gc.incid, ct.incid, hiv.prev, syph.incid, syphratio, syph.prev.hivpos, syph.prev.hivneg)

  return(out)
}


priors <- list(c("unif", 0.018, 0.030),
               c("unif", 2.0, 2.5),
               c("unif", 2.0, 2.5))


targets <- c(4.2, 6.6, 0.26, 0.9, 3.96, 0.103, 0.026)


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
