
library("EpiModelHIV")

sim <- NULL

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b4.rda scripts/burnin/")
system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataChosen.b4.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b4.rda")
dim(sim)

sim
cbind(sapply(sim, function(x) length(unique(x))))

# rgc, ugc, rct, uct
# targets <- c(0.083, 0.015, 0.118, 0.027)
targets <- c(0.102, 0.111, 0.141, 0.084)

nhist <- function(...) hist(..., col = "steelblue2", border = "white")

par(mfrow = c(2,2), mar = c(3,3,3,1), mgp = c(2,1,0))
nhist(sim$rgc.prev)
  abline(v = targets[1], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$rgc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$rgc.prev) - targets[1], 3)), cex = 0.8)
nhist(sim$ugc.prev)
  abline(v = targets[2], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$ugc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$ugc.prev) - targets[2], 3)), cex = 0.8)
nhist(sim$rct.prev)
  abline(v = targets[3], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$rct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$rct.prev) - targets[3], 3)), cex = 0.8)
nhist(sim$uct.prev)
  abline(v = targets[4], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(sim$uct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(sim$uct.prev) - targets[4], 3)), cex = 0.8)

rejection <- function(sim, targets, threshold) {
  diff.rgc <- abs(sim$rgc.prev - targets[1])
  diff.ugc <- abs(sim$ugc.prev - targets[2])
  diff.rct <- abs(sim$rct.prev - targets[3])
  diff.uct <- abs(sim$uct.prev - targets[4])

  choice <- which(diff.rgc <= threshold & diff.ugc <= threshold &
                    diff.rct <= threshold & diff.uct <= threshold)
  simChosen <- sim[choice, ]
  return(simChosen)
}

# round 1: 0.03
# round 2: 0.02
# round 3: 0.01
th <- 0.002

simChosen <- rejection(sim, targets = targets, threshold = th)

load("scripts/burnin/simDataChosen.b4.rda")
dim(simChosen)
simChosen

par(mfrow = c(2,2), mar = c(3,3,3,1), mgp = c(2,1,0))
  nhist(simChosen$rgc.prev)
  abline(v = targets[1], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$rgc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$rgc.prev) - targets[1], 3)), cex = 0.8)
nhist(simChosen$ugc.prev)
  abline(v = targets[2], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$ugc.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$ugc.prev) - targets[2], 3)), cex = 0.8)
nhist(simChosen$rct.prev)
  abline(v = targets[3], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$rct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$rct.prev) - targets[3], 3)), cex = 0.8)
nhist(simChosen$uct.prev)
  abline(v = targets[4], col = "firebrick", lty = 2, lwd = 2)
  abline(v = mean(simChosen$uct.prev), col = "darkblue", lwd = 2)
  mtext(paste0("diff=", round(mean(simChosen$uct.prev) - targets[4], 3)), cex = 0.8)


# for next batch

mn <- apply(simChosen, 2, mean)
sd <- apply(simChosen, 2, sd)
lo <- mn - sd
hi <- mn + sd

df <- data.frame(mn, sd, lo, hi)
round(df, 4)

save(simChosen, file = "scripts/burnin/simChosen.b3.rda")

file.rename("scripts/burnin/simDataAll.rda", "scripts/burnin/simDataAll.b3.rda")
system("scp scripts/burnin/simChosen.b3.rda hyak:/gscratch/csde/sjenness/sti/")

system("scp scripts/burnin/*.abcr.* hyak:/gscratch/csde/sjenness/sti/")
