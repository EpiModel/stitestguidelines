
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

                     rgc.tprob = x[2],
                     ugc.tprob = x[3],
                     rct.tprob = x[4],
                     uct.tprob = x[5],
                     syph.tprob = x[6],

                     # HIV acquisition
                     hiv.rgc.rr = x[7],
                     hiv.ugc.rr = x[8],
                     hiv.rct.rr = x[7],
                     hiv.uct.rr = x[8],
                     hiv.syph.rr = x[9],

                     syph.incub.sympt.prob = 0,
                     syph.prim.sympt.prob = 0.70,
                     syph.seco.sympt.prob = 0.85,
                     syph.earlat.sympt.prob = 0,
                     syph.latelat.sympt.prob = 0,
                     syph.tert.sympt.prob = 1.0,

                     syph.prim.sympt.prob.tx = 0.80,
                     syph.seco.sympt.prob.tx = 0.80,
                     syph.earlat.sympt.prob.tx = 0.10,
                     syph.latelat.sympt.prob.tx = 0.10,
                     syph.tert.sympt.prob.tx = 1.0,

                     ept.coverage = 0.0,
                     stianntest.gc.hivneg.coverage = 0.44,
                     stianntest.ct.hivneg.coverage = 0.44,
                     stianntest.syph.hivneg.coverage = 0.35,
                     stihighrisktest.gc.hivneg.coverage = 0.0,
                     stihighrisktest.ct.hivneg.coverage = 0.0,
                     stihighrisktest.syph.hivneg.coverage = 0.0,
                     stianntest.gc.hivpos.coverage = 0.61,
                     stianntest.ct.hivpos.coverage = 0.61,
                     stianntest.syph.hivpos.coverage = 0.57,
                     stihighrisktest.gc.hivpos.coverage = 0.0,
                     stihighrisktest.ct.hivpos.coverage = 0.0,
                     stihighrisktest.syph.hivpos.coverage = 0.0,

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

  gc.incid <- mean(tail(df$ir100.gc, 10))
  ct.incid <- mean(tail(df$ir100.ct, 10))
  hiv.prev <- mean(tail(df$i.prev, 10))
  syph.incid <- mean(tail(df$ir100.syph, 10))
  syph.prev <- mean(tail(df$prev.syph, 10))
  pssyph.prev <- mean(tail(df$prev.primsecosyph, 10))
  gctest.hivneg <- mean(tail(df$test.gc.12mo.hivneg, 10))
  gctest.hivpos <- mean(tail(df$test.gc.12mo.hivpos, 10))
  cttest.hivneg <- mean(tail(df$test.ct.12mo.hivneg, 10))
  cttest.hivpos <- mean(tail(df$test.ct.12mo.hivpos, 10))
  syphtest.hivneg <- mean(tail(df$test.syph.12mo.hivneg, 10))
  syphtest.hivpos <- mean(tail(df$test.syph.12mo.hivpos, 10))

  gcslope <- mean(df$ir100.gc[52] - df$ir100.gc[47])
  ctslope <- mean(df$ir100.ct[52] - df$ir100.ct[47])
  syphslope <- mean(df$ir100.syph[52] - df$ir100.syph[47])
  hivslope <- mean(df$ir100[52] - df$ir100[47])
  hivprevslope <- mean(df$i.prev[52] - df$i.prev[47])
  syphprevslope <- mean(df$prev.syph[52] - df$prev.syph[47])

  out <- c(gc.incid, ct.incid, hiv.prev, syph.incid,
           syph.prev, pssyph.prev,
           gctest.hivneg, gctest.hivpos,
           cttest.hivneg, cttest.hivpos,
           syphtest.hivneg, syphtest.hivpos,
           gcslope, ctslope, syphslope, hivslope,
           hivprevslope, syphprevslope)

  return(out)
}

priors <- list(c("unif", 0.50, 0.80), #rgc.tprob
               c("unif", 0.40, 0.60), #ugc.tprob
               c("unif", 0.27, 0.33), #rct.tprob
               c("unif", 0.20, 0.26), #uct.tprob
               c("unif", 0.01, 0.10), #syph.tprob
               c("unif", 1.77, 1.82), #rectal STI RR for HIV acquistion
               c("unif", 1.27, 1.34), #urethal STI RR for HIV acquistion
               c("unif", 1.60, 1.80)) #syph STI RR for HIV acquistion

# NG inc, CT inc, HIV prev, syph inc,
# syph prev, PS syph prev,
# NG test 12 months HIV neg, NG test 12 months HIV pos,
# CT test 12 months HIV neg, CT test 12 months HIV pos,
# Syph test 12 months HIV neg, Syph test 12 months HIV pos,
# GC inc slope, ct inc slope, syph inc slope, HIV inc slope,
# HIV prev slope, Syph prev slope
targets <- c(3.5, 5.6, 0.15, 2.6,
             0.012, 0.03,
             0.462, 0.641,
             0.458, 0.628,
             0.45, 0.68,
             0, 0, 0, 0,
             0, 0)

# NHBS NG/CT testing (Hoots 2011 and 2014 self-report data):
# 2014: NG: 46.2% HIV-MSM, 64.1%  HIV+ MSM
# 2014: CT: 45.8% HIV-MSM, 62.8%  HIV+ MSM
# NHBS syphilis testing (2014 self-report data Qian An): 45% HIV- MSM, 68% HIV+ MSM
# Flagg STD 2015 (MMP - 2008-2010 data among MSM): 54% syphilis, 20% NG, 20% CT
# Mattson CID 2017 (MMP- 2009-2013 data among MSM):
# Syphilis (2009-2013): 54, 57, 58, 60, 66
# NG/CT (2009-2013): 18, 22, 26, 32, 39

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
