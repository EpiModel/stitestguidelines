
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModelHIV"))

# Main Model Fx -----------------------------------------------------------

f <- function(x) {
  # set.seed(x[1])
  require(EpiModelHIV)
  load("est/nwstats.rda")
  load("est/fit.rda")

  param <- param_msm(nwstats = st,

                     rgc.tprob = 0.570, # x2
                     ugc.tprob = 0.500, # x3
                     rct.tprob = 0.245, # x4
                     uct.tprob = 0.205, # x5

                     rgc.asympt.int = 21*7, # x6
                     ugc.asympt.int = 21*7, # x6
                     rct.asympt.int = 45*7,
                     uct.asympt.int = 45*7,

                     rgc.sympt.prob = 0.16,
                     ugc.sympt.prob = 0.80,
                     rct.sympt.prob = 0.14,
                     uct.sympt.prob = 0.48,

                     hiv.rgc.rr = 1.97,
                     hiv.ugc.rr = 1.48,
                     hiv.rct.rr = 1.97,
                     hiv.uct.rr = 1.48,

                     stianntest.ct.hivneg.coverage = 0.44,
                     stianntest.ct.hivpos.coverage = 0.61,

                     stitest.start = 1,
                     sti.correlation.time = 0)

  init <- init_msm(nwstats = st,
                   prev.ugc = 0.0015,
                   prev.rgc = 0.0015,
                   prev.uct = 0.0015,
                   prev.rct = 0.0015,
                   prev.syph.B = 0,
                   prev.syph.W = 0)

  control <- control_msm(nsteps = 2600, verbose = FALSE)

  sim <- netsim(est, param, init, control)

  vars <- c("ir100.gc", "ir100.ct", "i.prev", "num")
  sim$epi <- sim$epi[which(names(sim$epi) %in% vars)]

  df <- as.data.frame(sim)

  last.year <- tail(df, 52)
  last.3years <- tail(df, 52*3)

  out <- c(mean(last.year$ir100.gc),
           mean(last.3years$ir100.gc),
           mean(last.year$ir100.ct),
           mean(last.3years$ir100.ct),
           mean(last.year$i.prev),
           mean(last.3years$i.prev))

}


# ABC Priors and Target Stats ---------------------------------------------

priors <- list(c("unif", 0.54, 0.60),
               c("unif", 0.47, 0.53),
               c("unif", 0.215, 0.275),
               c("unif", 0.175, 0.235),
               c("unif", 15*7, 25*7))

targets <- c(4.2, 4.2, 6.6, 6.6, 0.15, 0.15)


# Run ABC Prep ------------------------------------------------------------

prep <- abc_smc_prep(model = f,
                     prior = priors,
                     nb_simul = 100,
                     summary_stat_target = targets,
                     n_cluster = 28,
                     alpha = 0.5)
prep
saveRDS(prep, file = "data/abc.prep.rda")

# Batches for Wave 0
ceiling(prep$nb_simul/prep$n_cluster)

# Batches for Wave 1+
ceiling((prep$nb_simul - ceiling(prep$nb_simul * prep$alpha))/prep$n_cluster)

