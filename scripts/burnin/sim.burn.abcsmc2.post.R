
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

system("scp hyak:/gscratch/csde/sjenness/sti2/data/*.rda scripts/burnin/")

load("scripts/burnin/smc.fit.rda")
load("scripts/burnin/smc.fit.pacc2pct.250sim.rda")
load("scripts/burnin/smc.fit.pacc1pct.250sim.rda")


ls()
str(a)

p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

names(p) <- c("rgc.tprob", "ugc.tprob", "rct.tprob", "uct.tprob",
              "rgc.sympt.prob", "ugc.sympt.prob", "rct.sympt.prob", "uct.sympt.prob",
              "rgc.dur.asympt", "ugc.dur.asympt", "rct.dur.asympt", "uct.dur.asympt",
              "hiv.rect.rr", "hiv.ureth.rr")
names(s) <- c("rect.prev", "ureth.prev", "gc.incid", "ct.incid", "hiv.prev")

mean.s <- apply(s, 2, function(x) sum(x * w))
mean.p <- apply(p, 2, function(x) sum(x * w))

tar <- c(0.17, 0.07, 43, 48, 0.26)

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,2))
for (i in 1:ncol(s)) {
  hist(s[, i], col = "bisque2", border = "white", main = names(s)[i])
  abline(v = tar[i], lwd = 2, col = "red")
}

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(4,4))
for (i in 1:ncol(p)) {
  hist(p[, i], col = "bisque2", border = "white", main = names(p)[i])
}

save(mean.p, file = "scripts/burnin/abc.parms.v1.rda")

for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

