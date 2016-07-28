
## Analysis script for HIV risk compensation paper

library("EpiModelHIV")
library("dplyr")
source("scripts/analysis/fx.R")

# unlink("data/*")
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")

( fn <- list.files("data/", full.names = TRUE) )

for (i in seq_along(fn)) {
  load(fn[i])

  nsim <- sim$control$nsims
  cov <- rep(sim$param$prep.coverage, nsim)
  scint <- rep(sim$param$prep.sti.screen.int, nsim)
  rc <- rep(sim$param$rcomp.prob, nsim)

  epi <- truncate_sim(sim, at = 2600)$epi

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
  gc.prev <- as.numeric(tail(epi$prev.gc, 1))
  ct.prev <- as.numeric(tail(epi$prev.ct, 1))
  hiv.prev <- as.numeric(tail(epi$i.prev, 1))

  dft <- data.frame(cov, scint, rc, ir100.gc, ir100.ct, ir100.rgc, ir100.ugc,
                    ir100.rct, ir100.uct, tx.gc, tx.ct,
                    recov.rgc, recov.ugc, gc.prev, ct.prev, hiv.prev)
  if (i == 1) {
    df <- dft
  } else {
    df <- rbind(df, dft)
  }
  cat("*")
}

table(df$cov, df$scint, df$rc)

bycov <- group_by(df, scint, cov, rc)

summarise(bycov, mn = round(mean(ir100.gc), 1))
summarise(bycov, mn = round(mean(ir100.ct), 1))

summarise(bycov, mn = round(mean(ir100.rgc), 1))
summarise(bycov, mn = round(mean(ir100.ugc), 1))

summarise(bycov, mn = round(mean(hiv.prev), 3))

summarise(bycov, mn = round(mean(tx.gc), 1))
summarise(bycov, mn = round(mean(tx.ct), 1))

summarise(bycov, mn = round(mean(gc.prev), 3))
summarise(bycov, mn = round(mean(ct.prev), 3))

summarise(bycov, mn = round(mean(recov.rgc), 1))
summarise(bycov, mn = round(mean(recov.ugc), 1))


par(mfrow=c(1,1))
pal <- RColorBrewer::brewer.pal(3, "Set1")[1:2]
boxplot(ir100.gc ~ scint*cov, data = df, outline = FALSE, col = pal)




load(fn[2])
sim <- truncate_sim(sim, at = 2600)
plot(sim, y = "ir100.gc", mean.smooth = FALSE)

dat$epi$txGC <- rep(NA, length(dat$epi$num))
dat$epi$txCT <- rep(NA, length(dat$epi$num))



system("scp scripts/followup/*.fu.[Rs]* hyak:/gscratch/csde/sjenness/sti")
system("scp source/*.* hyak:/gscratch/csde/sjenness/sti/source")
