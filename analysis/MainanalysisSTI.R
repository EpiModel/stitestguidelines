
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
  dft <- data.frame(cov, scint, rc)

  epi <- truncate_sim(sim, at = 2600)$epi

  dft <- within(dft, {
    ir100.gc <- unname(colMeans(tail(epi$ir100.gc, 52)))
    ir100.ct <- unname(colMeans(tail(epi$ir100.ct, 52)))
    # ir100.rgc <- unname(colMeans(tail(epi$ir100.rgc, 52)))
    # ir100.ugc <- unname(colMeans(tail(epi$ir100.ugc, 52)))
    # ir100.rct <- unname(colMeans(tail(epi$ir100.rct, 52)))
    # ir100.uct <- unname(colMeans(tail(epi$ir100.uct, 52)))
    ir100 <- unname(colMeans(tail(epi$ir100, 52)))
    tx.gc <- unname(colMeans(tail(epi$txGC, 52)))
    tx.ct <- unname(colMeans(tail(epi$txCT, 52)))
    # recov.rgc <- unname(colMeans(tail(epi$recov.rgc, 52)))
    # recov.ugc <- unname(colMeans(tail(epi$recov.ugc, 52)))
    # recov.rct <- unname(colMeans(tail(epi$recov.rct, 52)))
    # recov.uct <- unname(colMeans(tail(epi$recov.uct, 52)))
    # rect.prev <- as.numeric(tail(epi$prev.rgcct, 1))
    # ureth.prev <- as.numeric(tail(epi$prev.ugcct, 1))
    # gc.prev <- as.numeric(tail(epi$prev.gc, 1))
    # ct.prev <- as.numeric(tail(epi$prev.ct, 1))
    hiv.prev <- as.numeric(tail(epi$i.prev, 1))
    # trans.main <- unname(colMeans(epi$trans.main, 52))
    # trans.casl <- unname(colMeans(epi$trans.casl, 52))
    # trans.inst <- unname(colMeans(epi$trans.inst, 52))
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


# Adding hazard ratios (denominator is the mean of baseline scenario)
bycov <- within(bycov, {
    haz.ir100 <- ir100 / colMeans(head(bycov[, "ir100"], 500))
    haz.ir100.gc <- ir100.gc / colMeans(head(bycov[, "ir100.gc"], 500))
    haz.ir100.ct <- ir100.ct / colMeans(head(bycov[, "ir100.ct"], 500))
})
#bycov <- group_by(df, asympt)

# Hazard ratios
summarise(bycov, hazard = round(mean(haz.ir100), 3),
          hazardql = round(quantile(haz.ir100, probs = 0.025), 3),
          hazardqu = round(quantile(haz.ir100, probs = 0.975), 3))
summarise(bycov, hazard = round(mean(haz.ir100.gc), 3),
          hazardql = round(quantile(haz.ir100.gc, probs = 0.025), 3),
          hazardqu = round(quantile(haz.ir100.gc, probs = 0.975), 3))
summarise(bycov, hazard = round(mean(haz.ir100.ct), 3),
          hazardql = round(quantile(haz.ir100.ct, probs = 0.025), 3),
          hazardqu = round(quantile(haz.ir100.ct, probs = 0.975), 3))

# Incidence, prevalence, and infections by partner type
summarise(bycov, mn = round(mean(hiv.prev), 3), ql = round(quantile(hiv.prev, probs = 0.025), 3) , qu = round(quantile(hiv.prev, probs = 0.975), 3))

summarise(bycov, mn = round(mean(ir100), 2), ql = round(quantile(ir100, probs = 0.025), 2) , qu = round(quantile(ir100, probs = 0.975), 2),
          transmainHIV = (sum(trans.main)/sum(incid)))

summarise(bycov, mn = round(mean(gc.prev), 3), ql = round(quantile(gc.prev, probs = 0.025), 3) , qu = round(quantile(gc.prev, probs = 0.975), 3))

summarise(bycov, mn = round(mean(ir100.gc), 2), ql = round(quantile(ir100.gc, probs = 0.025), 2) , qu = round(quantile(ir100.gc, probs = 0.975), 2),
          transGC = (sum(trans.main.gc)/sum(incid.gc)))

summarise(bycov, mn = round(mean(ct.prev), 3), ql = round(quantile(ct.prev, probs = 0.025), 3) , qu = round(quantile(ct.prev, probs = 0.975),3))

summarise(bycov, mn = round(mean(ir100.ct), 2), ql = round(quantile(ir100.ct, probs = 0.025), 2) , qu = round(quantile(ir100.ct, probs = 0.975), 2),
          transCT = (sum(trans.main.ct)/sum(incid.ct)))

summarise(bycov, mn = round(mean(sti.prev), 3), ql = round(quantile(sti.prev, probs = 0.025), 3) , qu = round(quantile(sti.prev, probs = 0.975), 3))

summarise(bycov, mn = round(mean(ir100.sti), 2), ql = round(quantile(ir100.sti, probs = 0.025), 2) , qu = round(quantile(ir100.sti, probs = 0.975), 2))


# Transmissions by partner type in last 500 observations - a different analysis
df <- tail(as.data.frame(sim), 500)
sum(df$trans.main) / sum(df$incid)
sum(df$trans.casl) / sum(df$incid)
sum(df$trans.inst) / sum(df$incid)

sum(df$trans.main.gc) / sum(df$incid.gc)
sum(df$trans.casl.gc) / sum(df$incid.gc)
sum(df$trans.inst.gc) / sum(df$incid.gc)

sum(df$trans.main.ct) / sum(df$incid.ct)
sum(df$trans.casl.ct) / sum(df$incid.ct)
sum(df$trans.inst.ct) / sum(df$incid.ct)


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
