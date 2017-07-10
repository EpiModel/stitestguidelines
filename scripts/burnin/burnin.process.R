
## Process burn-in
library("EpiModelHPC")
library("EpiModelHIV")

system("scp scripts/burnin/*.burn.[Rs]* hyak:/gscratch/csde/sjenness/sti")

# Examine output
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")

list.files("data/")

load("data/sim.n100.rda")

df <- as.data.frame(sim)

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", ylim = c(0.1, 0.3), qnts = 0.5, mean.lwd = 1)
abline(h = 0.26, lty = 2)
text(x = 0, y = 0.28, round(mean(tail(df$i.prev, 100)), 3))

plot(sim, y = "ir100.gc", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 4.2, lty = 2)
text(0, 1, round(mean(tail(df$ir100.gc, 520)), 2))

plot(sim, y = "ir100.ct", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 6.6, lty = 2)
text(0, 1, round(mean(tail(df$ir100.ct, 520)), 2))


sim <- truncate_sim(sim, at = 2600)
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", ylim = c(0.1, 0.3), qnts = 0.5, mean.lwd = 1)
abline(h = 0.26, lty = 2)
text(x = 0, y = 0.28, round(mean(tail(df$i.prev, 100)), 3))

plot(sim, y = "ir100.gc", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 4.2, lty = 2)
text(0, 1, round(mean(tail(df$ir100.gc, 52)), 2))

plot(sim, y = "ir100.ct", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 6.6, lty = 2)
text(0, 1, round(mean(tail(df$ir100.ct, 52)), 2))

save(sim, file = "data/sim.n100.rda")

# Other Calibration ---------------------------------------------------

# Merge sim files
sim <- merge_simfiles(simno = 120, indir = "data/", ftype = "max")

# Create function for selecting sim closest to target
mean_sim <- function(sim, targets) {

  nsims <- sim$control$nsims

  # Initialize distance vector
  dist <- rep(NA, nsims)

  # Obtain statistics and perform multivariable Euclidean distance calculation
  for (i in 1:nsims) {

      # Create data frame to draw statistics from
      df <- as.data.frame(x = sim, out = "vals", sim = i)

      # Create a vector of statistics
      calib <- c(mean(tail(df$ir100.gc, 10)),
                 mean(tail(df$ir100.ct, 10)),
                 mean(tail(df$ir100.syph, 10)),
                 mean(tail(df$i.prev, 10)),
                 mean(df$ir100.syph[5200] - df$ir100.syph[5190]),
                 mean(df$ir100.gc[5200] - df$ir100.gc[5190]),
                 mean(df$ir100.ct[5200] - df$ir100.ct[5190]))#,
                 #mean(df$ir100.syph[5190] - df$ir100.syph[5180]),
                 #mean(df$ir100.gc[5190] - df$ir100.gc[5180]),
                 #mean(df$ir100.ct[5190] - df$ir100.ct[5180]),
                 #mean(df$ir100.syph[5180] - df$ir100.syph[5170]),
                 #mean(df$ir100.gc[5180] - df$ir100.gc[5170]),
                 #mean(df$ir100.ct[5180] - df$ir100.ct[5170]))

      wts <- c(1, 1, 1, 1, 1, 1, 1)#, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

      # Iteratively calculate distance
      dist[i] <- sqrt(sum(((calib - targets)*wts)^2))
  }

  # Which sim minimizes distance
  meansim <- which.min(dist)
  return(meansim)
}

# Run function
mean_sim(sim, targets = c(3.5, 5.6, 2.6, 0.15, 0, 0, 0))#, 0, 0, 0, 0, 0, 0))


# Save burn-in file for FU sims
sim2 <- get_sims(sim, sims = 20)

par(mfrow = c(2,2), oma = c(0,0,2,0))
# plot(sim, y = "ir100")
# abline(h = 3.8, col = "red", lty = 2)
plot(sim, y = "i.prev", qnts = 0.90)
abline(h = 0.15, col = "red", lty = 2)
title("HIV Prevalence")
plot(sim, y = "ir100.gc", qnts = 0.90)
abline(h = 3.5, col = "red", lty = 2)
title("GC Incidence")
plot(sim, y = "ir100.ct", qnts = 0.90)
abline(h = 5.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph", qnts = 0.90)
abline(h = 2.6, col = "red", lty = 2)
title("Syph Incidence")
title("Summary of Sims", outer = TRUE)

tail(as.data.frame(sim2)$i.prev)
par(mfrow = c(2,2), oma = c(0,0,2,0))
# plot(sim, y = "ir100")
# abline(h = 3.8, col = "red", lty = 2)
plot(sim2, y = "i.prev")
abline(h = 0.15, col = "red", lty = 2)
title("HIV Prevalence")
plot(sim2, y = "ir100.gc")
abline(h = 3.5, col = "red", lty = 2)
title("GC Incidence")
plot(sim2, y = "ir100.ct")
abline(h = 5.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim2, y = "ir100.syph")
abline(h = 2.6, col = "red", lty = 2)
title("Syph Incidence")
title("Best-fitting Sim", outer = TRUE)

mean(tail(as.data.frame(sim2)$ir100.gc, 10))
mean(tail(as.data.frame(sim2)$ir100.ct, 10))
mean(tail(as.data.frame(sim2)$ir100.syph, 10))
mean(tail(as.data.frame(sim2)$i.prev, 10))

sim <- sim2
save(sim, file = "est/stimod.burnin.rda")
system("scp est/stimod.burnin.rda hyak:/gscratch/csde/sjenness/sti/est/")
