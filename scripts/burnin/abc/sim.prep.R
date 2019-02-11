
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModelHIV"))

# Main Model Fx -----------------------------------------------------------

f <- function(x) {
  set.seed(x[1])
  require(EpiModelHIV)
  load("est/nwstats.20k.rda")
  load("est/fit.20k.rda")

  param <- param_msm(nwstats = st,

                     rgc.tprob = x[2],
                     ugc.tprob = x[3],
                     rct.tprob = x[4],
                     uct.tprob = x[5],

                     rgc.asympt.rate = 1/(24.78753*7),
                     ugc.asympt.rate = 1/(24.78753*7),
                     rct.asympt.rate = 1/(44.28232*7),
                     uct.asympt.rate = 1/(44.28232*7),

                     rgc.sympt.prob = 0.16,
                     ugc.sympt.prob = 0.80,
                     rct.sympt.prob = 0.14,
                     uct.sympt.prob = 0.48,

                     hiv.rgc.rr = x[6], #1.97,
                     hiv.ugc.rr = x[7], #1.48,
                     hiv.rct.rr = x[6], #1.97,
                     hiv.uct.rr = x[7], #1.48,

                     ai.scale = x[8],
                     ai.scale.pospos = x[8],

                     stianntest.ct.hivneg.coverage = 0.44,
                     stianntest.ct.hivpos.coverage = 0.61,

                     stitest.start = 1,
                     sti.correlation.time = 0)

  init <- init_msm(nwstats = st,
                   prev.B = x[9], # 0.148,
                   prev.W = x[9], # 0.148,
                   prev.ugc = x[10], # 0.003,
                   prev.rgc = x[10], # 0.003,
                   prev.uct = x[11], # 0.005,
                   prev.rct = x[11], # 0.005,
                   prev.syph.B = 0,
                   prev.syph.W = 0)

  control <- control_msm(nsteps = 52*60,
                         verbose = FALSE)

  sim <- netsim(est, param, init, control)

  vars <- c("ir100.gc", "ir100.ct", "i.prev")
  sim$epi <- sim$epi[which(names(sim$epi) %in% vars)]

  years <- c(45:59)
  gc <- sapply(years, function(tt) mean(sim$epi$ir100.gc[(tt*52):(tt*52 + 51), ]))
  ct <- sapply(years, function(tt) mean(sim$epi$ir100.ct[(tt*52):(tt*52 + 51), ]))
  hiv <- sapply(years, function(tt) mean(sim$epi$i.prev[(tt*52):(tt*52 + 51), ]))

  out <- c(gc, ct, hiv)
  return(out)
}


# ABC Priors and Target Stats ---------------------------------------------

priors <- list(c("unif", 0.50, 0.57),
               c("unif", 0.40, 0.51),
               c("unif", 0.22, 0.27),
               c("unif", 0.18, 0.21),
               c("unif", 2.00, 2.30),
               c("unif", 1.45, 1.65),
               c("unif", 1.03, 1.08),
               c("unif", 0.146, 0.152),
               c("unif", 0.001, 0.005),
               c("unif", 0.003, 0.008))

years <- c(45:59)
targets <- c(4.2, 6.6, 0.15)
targets <- rep(targets, each = length(years))



# Run ABC Prep ------------------------------------------------------------

prep <- abc_smc_prep(model = f,
                     prior = priors,
                     nsims = 2500,
                     summary_stat_target = targets,
                     ncores = 28,
                     alpha = 0.2)
prep
saveRDS(prep, file = "scripts/burnin/abc/data/abc.prep.rda")

sbatch_master_abc(prep,
                  nwaves = 30,
                  master.file = "scripts/burnin/abc/master.sh",
                  mem = "100G",
                  ckpt = TRUE)
