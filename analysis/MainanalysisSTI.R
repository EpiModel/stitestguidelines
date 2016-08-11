
## Analysis script for HIV risk compensation paper

rm(list=ls())
library("EpiModelHIV")
library("dplyr")
source("analysis/fx.R")

unlink("data/*")
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")

( fn <- list.files("data/", full.names = TRUE) )

for (i in seq_along(fn)) {
  load(fn[i])

  nsim <- sim$control$nsims
  cov <- rep(sim$param$prep.coverage, nsim)
  scint <- rep(sim$param$prep.sti.screen.int, nsim)
  rc <- rep(sim$param$rcomp.prob, nsim)
  asympt <- rep(sim$param$gc.asympt.prob.tx, nsim)
  dft <- data.frame(cov, scint, rc, asympt)

  epi <- truncate_sim(sim, at = 2600)$epi

  dft <- within(dft, {
    ir100.gc <- unname(colMeans(tail(epi$ir100.gc, 52)))
    ir100.ct <- unname(colMeans(tail(epi$ir100.ct, 52)))
    ir100.rgc <- unname(colMeans(tail(epi$ir100.rgc, 52)))
    ir100.ugc <- unname(colMeans(tail(epi$ir100.ugc, 52)))
    ir100.rct <- unname(colMeans(tail(epi$ir100.rct, 52)))
    ir100.uct <- unname(colMeans(tail(epi$ir100.uct, 52)))
    tx.gc <- unname(colMeans(tail(epi$txGC, 52)))
    tx.ct <- unname(colMeans(tail(epi$txCT, 52)))
    recov.rgc <- unname(colMeans(tail(epi$recov.rgc, 52)))
    recov.ugc <- unname(colMeans(tail(epi$recov.ugc, 52)))
    recov.rct <- unname(colMeans(tail(epi$recov.rct, 52)))
    recov.uct <- unname(colMeans(tail(epi$recov.uct, 52)))
    rect.prev <- as.numeric(tail(epi$prev.rgcct, 1))
    ureth.prev <- as.numeric(tail(epi$prev.ugcct, 1))
    gc.prev <- as.numeric(tail(epi$prev.gc, 1))
    ct.prev <- as.numeric(tail(epi$prev.ct, 1))
    hiv.prev <- as.numeric(tail(epi$i.prev, 1))
  })

  if (i == 1) {
    df <- dft
  } else {
   df <- rbind(df, dft)
  }
  cat("*")
}

table(df$cov, df$asympt)

bycov <- group_by(df, scint, cov, rc)
bycov <- group_by(df, asympt)

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
