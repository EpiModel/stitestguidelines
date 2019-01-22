
library("EpiModelHIV")

list.files("data/")
# unlink(list.files("data/", full.names = TRUE))

load("data/sim.n100.rda")

print(sim)
# vars <- c("ir100.gc", "ir100.ct", "i.prev", "num")
# sim$epi <- sim$epi[which(names(sim$epi) %in% vars)]
df <- as.data.frame(sim)
names(df)

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "ir100.ct")
abline(h = 6.6)
plot(sim, y = "ir100.gc")
abline(h = 4.2)
plot(sim, y = "i.prev")
abline(h = 0.15)

plot(sim, y = "num")

# plot(sim, y = c("test.gc.12mo.hivdiag", "test.gc.12mo.nonhivdiag"))
# abline(h = c(0.458, 0.641), col = c("red", "blue"), lwd = 2, lty = 2)

# plot(sim, y = c("test.ct.12mo.hivdiag", "test.ct.12mo.nonhivdiag"))
# abline(h = c(0.458, 0.628), col = c("red", "blue"), lwd = 2, lty = 2)

