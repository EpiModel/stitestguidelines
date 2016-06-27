
library(methods)
suppressMessages(library(EpiModel))
library(doParallel)
library(foreach)

## Environmental Arguments
args <- commandArgs(trailingOnly = TRUE)
simno <- args[1]

load("est/fit.rda")

est <- est[[2]]

f <- function(est) {

  est$coef.form[1] <- est$coef.form[1] + runif(1, -0.15, 0.15)
  est$coef.form[2] <- est$coef.form[2] + runif(1, -0.15, 0.15)
  est$coef.form[3] <- est$coef.form[3] + runif(1, -0.15, 0.15)
  est$coef.form[4] <- est$coef.form[4] + runif(1, -0.15, 0.15)

  dx <- netdx(est, nsims = 1, nsteps = 300, verbose = FALSE,
              set.control.ergm = control.simulate.ergm(MCMC.burnin = 2e6))
  stat1 <- mean(tail(dx$stats[[1]][, "edges"], 100))
  stat2 <- mean(tail(dx$stats[[1]][, "concurrent"], 100))
  return(c(p1 = est$coef.form[1],
           p2 = est$coef.form[2],
           p3 = est$coef.form[3],
           p4 = est$coef.form[4],
           out1 = stat1,
           out2 = stat2))
}

registerDoParallel(16)
nsims <- 250
sout <- foreach(s = 1:nsims) %dopar% {
  f(est)
}

sim <- as.data.frame(do.call("rbind", sout))
save(sim, file = paste0("data/sim.", simno, ".rda"))

# rejection <- function(sim, target.stat = c(2022.5, 950), threshold = 0.05) {
#   diff1 <- abs(sim$out1 - target.stat[1])
#   diff2 <- abs(sim$out2 - target.stat[2])
#
#   diff <- diff1 * diff2
#   cutoff <- quantile(diff, threshold)
#
#   in.threshold <- which(diff <= cutoff)
#
#   post <- sim[in.threshold, ]
#   return(post)
# }
#
# post <- rejection(sim, threshold = 0.02)
# post
#
# colMeans(post)
#
# selection <- colMeans(post)
#
# est2 <- est
# est2$coef.form[1:4] <- selection[1:4]
#
# dx <- netdx(est2, nsteps = 300, nsims = 20, ncores = 20, dynamic = TRUE,
#             set.control.ergm = control.simulate.ergm(MCMC.burnin = 2e6))
# dx2
# plot(dx)

