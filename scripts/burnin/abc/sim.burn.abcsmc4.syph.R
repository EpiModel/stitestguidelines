
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
                     ai.scale = 1.03,

                     syph.earlat.rr = 0.5,
                     incu.syph.int = 27,
                     prim.syph.int = 60,
                     seco.syph.int = 120,
                     earlat.syph.int = 365 - 27 - 60 - 120,
                     latelat.syph.int = 9 * 52 * 7,
                     latelatelat.syph.int = 20 * 52 * 7,
                     tert.syph.int = 20 * 52 * 7,
                     syph.tert.prog.prob = 0.00015625599,

                     rgc.tprob = x[2],
                     ugc.tprob = x[3],
                     rct.tprob = x[4],
                     uct.tprob = x[5],
                     syph.tprob = x[6],


                     # HIV acquisition
                     hiv.rgc.rr = 1.80292790,
                     hiv.ugc.rr = 1.1989083,
                     hiv.rct.rr = 1.80292790,
                     hiv.uct.rr = 1.1989083,
                     hiv.syph.rr = 1.62,

                     # HIV transmission
                     hiv.trans.gc.rr = 1,
                     hiv.trans.ct.rr = 1,
                     hiv.trans.syph.rr = 1,

                     syph.prim.sympt.prob.tx = 0.60,
                     syph.seco.sympt.prob.tx = 0.688235,
                     syph.earlat.sympt.prob.tx = 0.10,
                     syph.latelat.sympt.prob.tx = 0.10,
                     syph.tert.sympt.prob.tx = 1.0,

                     syph.prim.asympt.prob.tx = 1,
                     syph.seco.asympt.prob.tx = 1,
                     syph.earlat.asympt.prob.tx = 1,
                     syph.latelat.asympt.prob.tx = 1,
                     syph.tert.asympt.prob.tx = 1,

                     hivdx.syph.sympt.tx.rr = 1.5,

                     prep.coverage = 0.0,
                     ept.coverage = 0.0,
                     stianntest.coverage = 0.1,
                     stihighrisktest.coverage = 0.1,

                     prep.start = 7000,
                     stitest.start = 5201,
                     ept.start = 7000,

                     stitest.elig.model = "sti",

                     stitest.active.int = 364,
                     sti.highrisktest.int = 182,
                     ept.risk.int = 60)

  init <- init_msm(nwstats = st)

  control <- control_msm(simno = 1,
                         nsteps = 5200,
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


priors <- list(c("unif", 0.44, 0.48),
               c("unif", 0.33, 0.36),
               c("unif", 0.20, 0.22),
               c("unif", 0.19, 0.20),
               c("unif", 0.145, 0.147))

targets <- c(3.5, 5.6, 0.15, 2.6, 0, 0, 0, 0, 0)

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
