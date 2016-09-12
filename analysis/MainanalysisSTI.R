
## Analysis script for HIV risk compensation paper

rm(list=ls())
library("EpiModelHIV")
library("dplyr")
source("analysis/fx.R")

unlink("data/*")
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")

( fn <- list.files("data/", full.names = TRUE) )

# CROI Table ----------------------------------------------------------

sim.base <- merge_simfiles(1000)
sim.comp <- merge_simfiles(1001)

epi_stats(sim.base, at = 520, qnt.low = 0.025, qnt.high = 0.975)
epi_stats(sim.base, sim.comp, at = 520, qnt.low = 0.025, qnt.high = 0.975)





for (i in seq_along(fn)) {
  load(fn[i])

  nsim <- sim$control$nsims
  cov <- rep(sim$param$prep.coverage, nsim)
  scint <- rep(sim$param$prep.sti.screen.int, nsim)
  rc <- rep(sim$param$rcomp.prob, nsim)
  dft <- data.frame(cov, scint, rc)

  epi <- truncate_sim(sim, at = 2600)$epi

  dft <- within(dft, {
    ir100.gc <- unname(colMeans(tail(epi$ir100.gc, 52)))
    ir100.ct <- unname(colMeans(tail(epi$ir100.ct, 52)))
    ir100 <- unname(colMeans(tail(epi$ir100, 52)))
    tx.gc <- unname(colMeans(tail(epi$txGC, 52)))
    tx.ct <- unname(colMeans(tail(epi$txCT, 52)))
    hiv.prev <- as.numeric(tail(epi$i.prev, 1))
  })

  if (i == 1) {
    df <- dft
  } else {
   df <- rbind(df, dft)
  }
  cat("*")
}

table(df$cov, df$scint)

bycov <- group_by(df, scint)

summarise(bycov, mn = round(mean(hiv.prev), 3))
summarise(bycov, mn = round(mean(ir100), 2))
summarise(bycov, mn = round(mean(ir100.gc), 2))
summarise(bycov, mn = round(mean(ir100.ct), 2))

par(mar = c(3,3,1,1), mgp = c(2,1,0))
sim <- truncate_sim(sim, at = 2600)
plot(sim, y = "i.prev", ylim = c(0.1, 0.3), qnts = 0.5, mean.lwd = 1)
plot(sim, y = "ir100", ylim = c(0, 6), mean.smooth = TRUE, mean.lwd = 1)
plot(sim, y = "ir100.gc", mean.smooth = TRUE, mean.lwd = 1, ylim = c(0, 10))
plot(sim, y = "ir100.ct", mean.smooth = TRUE, mean.lwd = 1, ylim = c(0, 10))




# Other summary statistics
summarise(bycov, mn = round(mean(hiv.prev), 3))
summarise(bycov, mn = round(mean(gc.prev), 3))
summarise(bycov, mn = round(mean(ct.prev), 3))
summarise(bycov, mn = round(mean(rect.prev), 3))
summarise(bycov, mn = round(mean(ureth.prev), 3))

summarise(bycov, mn = round(mean(tx.gc), 1))
summarise(bycov, mn = round(mean(tx.ct), 1))

summarise(bycov, mn = round(mean(ir100.gc), 1))
summarise(bycov, mn = round(mean(ir100.ct), 1))

summarise(bycov, mn = round(mean(ir100.rgc), 1))
summarise(bycov, mn = round(mean(ir100.ugc), 1))

summarise(bycov, mn = round(mean(ir100.rct), 1))
summarise(bycov, mn = round(mean(ir100.uct), 1))

summarise(bycov, mn = round(mean(recov.rgc), 1))
summarise(bycov, mn = round(mean(recov.ugc), 1))
summarise(bycov, mn = round(mean(recov.rct), 1))
summarise(bycov, mn = round(mean(recov.uct), 1))

names(df)
names(sim$epi)

par(mfrow=c(1,1), mar = c(3,3,1,1), mgp = c(2,1,0))
boxplot(ir100.gc ~ cov, data = df, outline = FALSE, col = "royalblue1")
boxplot(ureth.prev ~ cov, data = df, outline = FALSE, col = "royalblue1")

boxplot(ir100.gc ~ asympt, data = df, outline = FALSE, col = "royalblue1")
boxplot(ir100.ct ~ asympt, data = df, outline = FALSE, col = "royalblue1")

boxplot(gc.prev ~ asympt, data = df, outline = FALSE, col = "royalblue1")
boxplot(ct.prev ~ asympt, data = df, outline = FALSE, col = "royalblue1")


pal <- brewer_ramp(length(fn), "Blues")
for (i in seq_along(fn)) {
  load(fn[i])
  sim <- truncate_sim(sim, at = 2600)
  plot(sim, y = "prev.ct", mean.col = pal[i], qnts = FALSE, qnts.col = pal[i],
       add = i > 1)
}
