
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

plot(sim, y = "ir100", ylim = c(0, 6), mean.smooth = TRUE, mean.lwd = 1)
abline(h = 3.8, lty = 2)
text(0, 1, round(mean(tail(df$ir100, 260)), 2))

plot(sim, y = "ir100.gc", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 4.2, lty = 2)
text(0, 1, round(mean(tail(df$ir100.gc, 260)), 2))

plot(sim, y = "ir100.ct", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 6.6, lty = 2)
text(0, 1, round(mean(tail(df$ir100.ct, 260)), 2))





# Other Calibration ---------------------------------------------------

# Set targets for multivariate Euclidean distance calculation
# rect.prev, ureth.prev, gc.incid, ct.incid, hiv.incid, hiv.prev
targets <- c(0.135, 0.046, 23.2, 26.8, 3.8, 0.26)

# Merge sim files
sim <- merge_simfiles(simno = 100, indir = "data/", ftype = "max")

# Create function for selecting sim closest to target
mean_sim <- function(sim, targets) {

  nsims <- sim$control$nsims

  # Initialize distance vector
  dist <- c()
  calib <- c()
  targets <- as.data.frame(targets)

  # Obtain statistics and perform multivariable Euclidean distance calculation in a for-loop
  for (i in 1:nsims) {

      # Create data frame to draw statistics from
      df <- as.data.frame(x = sim, out = "vals", sim = i)

      # Create a vector of statistics
      calib[i] <- as.data.frame(c(mean(tail(df$ir100.gc, 52)),
                                  mean(tail(df$ir100.ct, 52)),
                                  mean(tail(df$i.prev, 52))))

      # Iteratively calculate distance
      dist[i] <- sqrt(sum((calib[i] - targets)^2))
  }

  # Which sim minimizes distance
  meansim <- which.min(dist)
  return(meansim)
}

# Run function
mean_sim(sim, targets = c(4.2, 6.6, 0.26))


# Save burn-in file for FU sims
sim <- get_sims(sim, sims = 115)
tail(as.data.frame(sim)$i.prev)
mean(tail(as.data.frame(sim)$ir100.gc, 52))
mean(tail(as.data.frame(sim)$ir100.ct, 52))

save(sim, file = "est/stimod.burnin.rda")
system("scp hyak:/gscratch/csde/sjenness/sti/est/stimod.burnin.rda est/")
