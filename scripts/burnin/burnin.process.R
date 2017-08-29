
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
sim <- merge_simfiles(simno = 109, indir = "data/", ftype = "max")

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
                 mean(df$ir100.ct[5200] - df$ir100.ct[5190]),
                 mean(df$ir100.syph[5190] - df$ir100.syph[5180]),
                 mean(df$ir100.gc[5190] - df$ir100.gc[5180]),
                 mean(df$ir100.ct[5190] - df$ir100.ct[5180]))
                 # mean(df$ir100.syph[5180] - df$ir100.syph[5170]),
                 # mean(df$ir100.gc[5180] - df$ir100.gc[5170]),
                 # mean(df$ir100.ct[5180] - df$ir100.ct[5170]))

      wts <- c(4, 4, 4, 4, 3, 3, 3, 1, 1, 1)#, 1, 1, 1)

      # Iteratively calculate distance
      dist[i] <- sqrt(sum(((calib - targets)*wts)^2))
  }

  # Which sim minimizes distance
  meansim <- which.min(dist)
  return(meansim)
}

# Run function
mean_sim(sim, targets = c(3.5, 5.6, 2.6, 0.15, 0, 0, 0, 0, 0, 0))#, 0, 0, 0, 0, 0, 0))


# Save burn-in file for FU sims
sim2 <- get_sims(sim, sims = 14)


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
# plot(sim2, y = "ir100")
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

par(mfrow = c(2,2), oma = c(0,0,2,0))
#plot(sim, y = "stiasympttests", mean.col = "blue", qnts.col = "blue", qnts = 0.5, ylab = "# of Tests")
plot(sim, y = "GCasympttests.hivneg", mean.col = "red", qnts.col = "red", qnts = 0.5)
plot(sim, y = "GCasympttests.hivpos", mean.col = "blue", qnts.col = "blue", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1), col = c("red", "blue"),
       legend = c("NG-non-HIV dx", "NG-HIV dx"))
plot(sim, y = "CTasympttests.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5)
plot(sim, y = "CTasympttests.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1), col = c("green", "orange"),
       legend = c("CT-non-HIV dx", "CT-HIV dx"))
plot(sim, y = "syphasympttests.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5)
plot(sim, y = "syphasympttests.hivpos", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1), col = c("purple", "brown"),
       legend = c("Syph-non-HIV dx", "Syph-HIV dx"))
#title("STI Testing - 10% Coverage \n Serostatus-specific" )
title("STI Testing - Serostatus-Specific", outer = TRUE)

plot(sim2, y = "stiasympttests", mean.col = "blue", ylab = "# of Tests")
plot(sim2, y = "GCasympttests", mean.col = "red", add = TRUE)
plot(sim2, y = "CTasympttests", mean.col = "green", add = TRUE)
plot(sim2, y = "syphasympttests", mean.col = "purple", add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1), col = c("blue", "red", "green", "purple"),
       legend = c("STI", "NG", "CT", "Syph"))
title("STI Testing", outer = TRUE)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph", qnts = 0.90)
#abline(h = 0.01, lty = c(2), col = 'red')
title("P&S Syphilis Prevalence")
plot(sim, y = "prev.syph", qnts = 0.90)
title("Syphilis (All Stages) Prevalence")
plot(sim, y = "early.late.syphratio", ylim = c(0, 1.0))
title("Ratio of Early to Late \n Syphilis Cases")
#abline(h = 0.2, lty = c(2), col = 'red')
plot(sim, y = "early.late.diagsyphratio", ylim = c(0, 1.0))
title("Ratio of Diagnosed Early to Late \n Syphilis Cases")
#abline(h = 0.5, lty = c(2), col = 'red')
title("Syphilis Prevalence Measures", outer = TRUE)

par(mfrow = c(1, 2))
plot(sim, y = "prev.stage.incub", mean.col = "green", qnts.col = "green",
     qnts = 0.5, qnts.alpha = 0.5, ylim = c(0, 0.4), ylab = "Proportion of all prevalent infections")
plot(sim, y = "prev.stage.prim", mean.col = "blue", qnts.col = "blue",
     qnts = 0.5, qnts.alpha = 0.5, add = TRUE)
plot(sim, y = "prev.stage.seco", mean.col = "red", qnts.col = "red",
     qnts = 0.5, qnts.alpha = 0.5, add = TRUE)
plot(sim, y = "prev.stage.earlat", mean.col = "orange", qnts.col = "orange",
     qnts = 0.5, qnts.alpha = 0.5, add = TRUE)
legend("topleft", lty = c(1,1,1,1), col = c("green", "blue", "red", "orange"),
       legend = c("Incub", "Prim", "Seco", "Early Latent"))

plot(sim, y = "prev.stage.latelat", mean.col = "purple", qnts.col = "purple",
     qnts = 0.5, ylab = "Proportion of all prevalent infections")
plot(sim, y = "prev.stage.tert", mean.col = "black", qnts.col = "black",
     qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1,1), col = c("purple", "black"),
       legend = c("Late Latent", "Tertiary"))
title("Prevalence of Stage-Specific Syphilis", outer = TRUE)

df <- as.data.frame(x = sim, out = "vals")
sum(tail(df$num.newearlydiagsyph, 52))
sum(tail(df$num.newlatediagsyph, 52))

mean(tail(as.data.frame(sim2)$ir100.gc, 10))
mean(tail(as.data.frame(sim2)$ir100.ct, 10))
mean(tail(as.data.frame(sim2)$ir100.syph, 10))
mean(tail(as.data.frame(sim2)$i.prev, 10))
mean(tail(as.data.frame(sim2)$ir100.gc, 5))
mean(tail(as.data.frame(sim2)$ir100.ct, 5))
mean(tail(as.data.frame(sim2)$ir100.syph, 5))
mean(tail(as.data.frame(sim2)$i.prev, 5))

sim <- sim2
save(sim, file = "est/stimod.burnin.rda")
system("scp est/stimod.burnin.rda hyak:/gscratch/csde/sjenness/sti/est/")
