
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
                     syph.prim.sympt.prob = x[10],
                     syph.seco.sympt.prob = x[11],
                     syph.earlat.sympt.prob = 0,
                     syph.latelat.sympt.prob = 0,
                     syph.tert.sympt.prob = 1.0,

                     syph.prim.sympt.prob.tx = x[12],
                     syph.seco.sympt.prob.tx = x[13],
                     syph.earlat.sympt.prob.tx = 0.10,
                     syph.latelat.sympt.prob.tx = 0.10,
                     syph.tert.sympt.prob.tx = 1.0,

                     ept.coverage = 0.0,
                     stianntest.gc.hivneg.coverage = x[14],
                     stianntest.ct.hivneg.coverage = x[15],
                     stianntest.syph.hivneg.coverage = x[16],
                     stihighrisktest.gc.hivneg.coverage = 0.0,
                     stihighrisktest.ct.hivneg.coverage = 0.0,
                     stihighrisktest.syph.hivneg.coverage = 0.0,
                     stianntest.gc.hivpos.coverage = x[17],
                     stianntest.ct.hivpos.coverage = x[18],
                     stianntest.syph.hivpos.coverage = x[19],
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

  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.prev <- mean(df$i.prev)
  syph.incid <- mean(df$ir100.syph)
  syph.prev <- mean(df$prev.syph)
  pssyph.prev <- mean(df$prev.primsecosyph)
  gctest.hivneg <- mean(df$test.gc.12mo.hivneg)
  gctest.hivpos <- mean(df$test.gc.12mo.hivpos)
  cttest.hivneg <- mean(df$test.ct.12mo.hivneg)
  cttest.hivpos <- mean(df$test.ct.12mo.hivpos)
  syphtest.hivneg <- mean(df$test.syph.12mo.hivneg)
  syphtest.hivpos <- mean(df$test.syph.12mo.hivpos)

  gcslope <- mean(df$ir100.gc[52] - df$ir100.gc[47])
  ctslope <- mean(df$ir100.ct[52] - df$ir100.ct[47])
  syphslope <- mean(df$ir100.syph[52] - df$ir100.syph[47])
  hivslope <- mean(df$ir100[52] - df$ir100[47])
  hivprevslope <- mean(df$i.prev[52] - df$i.prev[47])
  syphprevslope <- mean(df$prev.syph[52] - df$prev.syph[47])

  out <- c(gc.incid, ct.incid, hiv.prev, syph.incid,
           syph.prev, pssyph.prev,
           gctest.hivneg, gctest.hivpos, cttest.hivneg,
           cttest.hivpos, syphtest.hivneg, syphtest.hivpos,
           gcslope, ctslope, syphslope, hivslope,
           hivprevslope, syphprevslope)

  return(out)
}


priors <- list(c("unif", 0.40, 0.48), #rgc.tprob
               c("unif", 0.33, 0.40), #ugc.tprob
               c("unif", 0.22, 0.28), #rct.tprob
               c("unif", 0.17, 0.22), #uct.tprob
               c("unif", 0.20, 0.24), #syph.tprob
               c("unif", 1.70, 2.00), #rectal STI RR for HIV acquistion
               c("unif", 1.20, 1.40), #urethal STI RR for HIV acquistion
               c("unif", 1.50, 2.00), #syph STI RR for HIV acquistion
               c("unif", 0.60, 0.80), #syph.prob.sympt.prob
               c("unif", 0.70, 0.90), #syph.seco.sympt.prob
               c("unif", 0.75, 0.90), #syph.prim.sympt.prob.tx
               c("unif", 0.75, 0.90), #syph.seco.sympt.prob.tx
               c("unif", 0.40, 0.45), #stianntest.gc.hivneg.coverage
               c("unif", 0.40, 0.45), #stianntest.ct.hivneg.coverage
               c("unif", 0.10, 0.45), #stianntest.syph.hivneg.coverage
               c("unif", 0.60, 0.65), #stianntest.gc.hivpos.coverage
               c("unif", 0.60, 0.65), #stianntest.ct.hivpos.coverage
               c("unif", 0.65, 0.70)) #stianntest.syph.hivpos.coverage

# NG inc, CT inc, HIV prev, syph inc, syph prev, PS syph prev,
# NG test 12 months HIV neg, NG test 12 months HIV pos,
# CT test 12 months HIV neg, CT test 12 months HIV pos,
# Syph test 12 months HIV neg, Syph test 12 months HIV pos,
# GC inc slope, ct inc slope, syph inc slope, HIV inc slope,
# HIV prev slope, Syph prev slope
targets <- c(3.5, 5.6, 0.15, 2.6, 0.02, 0.01, 0.462, 0.641, 0.458, 0.628,
             0.45, 0.68, 0, 0, 0, 0, 0, 0)
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
