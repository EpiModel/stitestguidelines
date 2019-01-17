
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

# Main Model Fx -----------------------------------------------------------

myfunc <- function(x) {
  set.seed(x[1])
  require(EpiModel)
  load("fit.rda")
  est <- est[[2]]
  est$coef.form[1:4] <- est$coef.form[1:4] + x[2:5]
  dx <- netdx(est, nsims = 1, nsteps = 500, verbose = FALSE)
  stats <- get_nwstats(dx)[, 3:6]
  out <- unname(colMeans(tail(stats, 100)))
  return(out)
}


# ABC Priors and Target Stats ---------------------------------------------

priors <- list(c("unif", -1, 0.2),
               c("unif", -0.5, 0.5),
               c("unif", -0.5, 0.5),
               c("unif", -0.5, 0.5))

load("fit.rda")
est <- est[[2]]
targets <- est$target.stats[1:4]



# Run ABC Prep ------------------------------------------------------------

prep <- abc_smc_prep(model = myfunc,
                     prior = priors,
                     nb_simul = 1000,
                     summary_stat_target = targets,
                     n_cluster = 16,
                     alpha = 0.2)
prep
saveRDS(prep, file = "data/abc.prep.rda")

# Batches for Wave 0
ceiling(prep$nb_simul/prep$n_cluster)

# Batches for Wave 1+
ceiling((prep$nb_simul - ceiling(prep$nb_simul * prep$alpha))/prep$n_cluster)
