
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

names(p) <- c("rsyph.tprob", "usyph.tprob", "hiv.rsyph.rr", "hiv.usyph.rr", "rgc.tprob","ugc.tprob",
              "rct.tprob","uct.tprob", "hiv.rsti.rr", "hiv.usti.rr")

names(s) <- c("gc.incid", "ct.incid", "hiv.incid", "syph.incid","hiv.prev",
              "gcslope", "ctslope", "syphslope", "hivslope")

comb <- cbind(s, p)
View(comb)

( mean.s <- apply(s, 2, function(x) sum(x * w)) )
( mean.p <- apply(p, 2, function(x) sum(x * w)) )

# hist(s$prev.primsecosyph.hivpos / s$prev.primsecosyph.hivneg)

tar <- c(4.2, 6.6, 3.8, 2.0, 0.26)

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,2))
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
