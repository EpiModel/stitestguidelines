
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

system("scp hyak:/gscratch/csde/sjenness/sti2/data/*.rda scripts/burnin/")

## this was the first batch of fits to PrEP demo project

# 100 sim at 5%
load("scripts/burnin/smc.fit.rda")

# 250 sim at 2%
load("scripts/burnin/smc.fit.pacc2pct.250sim.rda")

# 250 sim at 1%
load("scripts/burnin/smc.fit.pacc1pct.250sim.rda")

## second batch of fits to PrEP demo project
##    main correction was getting STIs started at 5% to fix GC burnin issues
load("scripts/burnin/smc.2pct.250sim.rda")


## batch of fits to ATL::Involvement data
load("scripts/burnin/smc.atl.raceavg.5pct.100sim.rda")


## averaged ATL/demo fits
load("scripts/burnin/abc/smc.avg.5pct.100sim.rda")

p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

names(p) <- c("rgc.tprob", "ugc.tprob", "rct.tprob", "uct.tprob",
              "rgc.sympt.prob", "ugc.sympt.prob", "rct.sympt.prob", "uct.sympt.prob",
              "rgc.dur.asympt", "ugc.dur.asympt", "rct.dur.asympt", "uct.dur.asympt",
              "hiv.rect.rr", "hiv.ureth.rr", "prob.cease")

# for PrEP demo project fits
names(s) <- c("rect.prev", "ureth.prev", "gc.incid", "ct.incid", "hiv.prev")

# for ATL fits
names(s) <- c("rgc.prev", "ugc.prev", "rct.prev", "uct.prev",
              "rgc.incid", "ugc.incid", "rct.incid", "uct.incid",
              "hiv.prev")

# for averaged fits
names(s) <- c("rect.prev", "ureth.prev", "gc.incid", "ct.incid", "hiv.incid", "hiv.prev")


( mean.s <- apply(s, 2, function(x) sum(x * w)) )
( mean.p <- apply(p, 2, function(x) sum(x * w)) )


tar.demo <- c(0.17, 0.07, 43, 48, 0.26)
tar.atl <- c(0.083, 0.015, 0.118, 0.027,
             6.19, 1.07, 7.81, 3.75,
             0.26)
tar.avg <- c(0.135, 0.046, 23.2, 26.8, 3.8, 0.26)

data.frame(mean.s, tar.atl)

mean.p

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,2))
for (i in 1:ncol(s)) {
  hist(s[, i], col = "bisque2", border = "white", main = names(s)[i])
  abline(v = tar[i], lwd = 2, col = "red")
}

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(4,4))
for (i in 1:ncol(p)) {
  hist(p[, i], col = "bisque2", border = "white", main = names(p)[i])
}

save(mean.p, file = "scripts/burnin/abc/abc.avg.parms.5pct.rda")

for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

