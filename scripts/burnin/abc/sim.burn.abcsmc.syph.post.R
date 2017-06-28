
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

names(p) <- c("rsyph.tprob", "usyph.tprob", "rgc.tprob","ugc.tprob",
              "rct.tprob","uct.tprob", "hiv.rsti.rr", "hiv.usti.rr")

names(s) <- c("gc.incid", "ct.incid", "hiv.prev","syph.incid",
              "gcslope", "ctslope", "syphslope", "hivslope", "hivprevslope")

comb <- cbind(s, p)
#View(comb)

( mean.s <- apply(s, 2, function(x) sum(x * w)) )
( mean.p <- apply(p, 2, function(x) sum(x * w)) )

# hist(s$prev.primsecosyph.hivpos / s$prev.primsecosyph.hivneg)

tar <- c(3.5, 5.0, 2.0, 0.15, 0, 0, 0, 0, 0)

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,3))
for (i in 1:ncol(s)) {
  hist(s[, i], col = "bisque2", border = "white", main = names(s)[i])
  abline(v = tar[i], lwd = 2, col = "red")
}

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,4))
for (i in 1:ncol(p)) {
  hist(p[, i], col = "bisque2", border = "white", main = names(p)[i])
  abline(v = mean.p[i], lwd = 2, col = "red")
}

for (i in seq_along(mean.p)) {
    assign(names(mean.p)[i], unname(mean.p[i]))
}

save(mean.p, file = "scripts/burnin/abc/abc.avg.parms.syph.100sim.1pct.rda")
save(mean.p, file = "est/abc.syph.parms.rda")
