
library("EpiModelHIV")

sim <- NULL

# All data ------------------------------------------------------------

system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataAll.b1.rda scripts/burnin/")
load("scripts/burnin/simDataAll.b1.rda")
dim(sim)

# sim
cbind(sapply(sim, function(x) length(unique(x))))

# rgc, ugc, rct, uct
# targets <- c(0.083, 0.015, 0.118, 0.027)
# targets <- c(0.102, 0.111, 0.141, 0.084)

# rect.prev, ureth.prev, gc.incid, ct.incid, hiv.prev
targets <- c(0.17, 0.07, 43, 48, 0.26)

nhist <- function(...) hist(..., col = "seagreen3", border = "white")

par(mfrow = c(2,3), mar = c(3,3,3,1), mgp = c(2,1,0))
stats <- 15:19
for (i in seq_along(stats)) {
  nhist(sim[[stats[i]]], main = names(sim[stats[i]]))
  abline(v = mean(sim[[stats[i]]]), col = "blue", lwd = 2)
  abline(v = median(sim[[stats[i]]]), col = "blue", lty = 2, lwd = 2)
  abline(v = targets[i], col = "red", lwd = 2)
}

par(mfrow = c(4,4))
for (i in 1:14) {
  nhist(sim[[i]], main = names(sim[i]))
}


# Accepted data -------------------------------------------------------



# batch 1: 0.005
threshold <- 0.005

simChosen <- rejection(sim, targets = targets, threshold = 0.015)
simChosen <- rejection2(sim, targets = targets, threshold = 0.01)

# system("scp hyak:/gscratch/csde/sjenness/sti/data/simDataChosen.b4.rda scripts/burnin/")
# load("scripts/burnin/simDataChosen.b4.rda")
# simChosen

par(mfrow = c(2,3), mar = c(3,3,3,1), mgp = c(2,1,0))
stats <- 15:19
for (i in seq_along(stats)) {
  nhist(simChosen[[stats[i]]], main = names(simChosen[stats[i]]))
  abline(v = mean(simChosen[[stats[i]]]), col = "blue", lwd = 2)
  abline(v = median(simChosen[[stats[i]]]), col = "blue", lty = 2, lwd = 2)
  abline(v = targets[i], col = "red", lwd = 2)
}

mn <- apply(sim, 2, mean)
sds <- apply(sim, 2, sd)
lo <- mn - sds
hi <- mn + sds
data.frame(mn, lo, hi)

mn <- apply(simChosen, 2, mean)
sds <- apply(simChosen, 2, sd)
lo <- mn - sds
hi <- mn + sds
data.frame(mn, lo, hi)

par(mfrow = c(3,4))
for (i in 1:12) {
  nhist(simChosen[[i]], main = names(simChosen[i]))
}
for (i in 1:12) {
  plot(density(simChosen[[i]]), main = names(simChosen[i]))
}

for (i in 1:12) {
  d <- density(simChosen[[i]])
  cat(names(simChosen[i]), ":", d$x[which.max(d$y)], "\n")
}
pairs(simChosen[,1:12])

# for next batch
save(simChosen, file = "scripts/burnin/simChosen.b1.rda")

system("scp scripts/burnin/simChosen.b1.rda hyak:/gscratch/csde/sjenness/sti/")

system("scp scripts/burnin/*.abcr.* hyak:/gscratch/csde/sjenness/sti/")
system("scp source/*.* hyak:/gscratch/csde/sjenness/sti/source/")
