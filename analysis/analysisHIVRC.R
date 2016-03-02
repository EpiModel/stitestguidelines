
## Analysis script for HIV risk compensation paper

library(EpiModelHIVmsm)
source("scripts/analysis/fx.R")

# scp hyak:/gscratch/csde/camp/data/save/*.rda data/

steps <- 52 *  10

## Base model: sim.n2000
load("data/sim.n2000.rda")
sim.base <- truncate_sim(sim, at = 2600)

mn <- as.data.frame(sim.base)
ql <- as.data.frame(sim.base, out = "qnt", qval = 0.05)
qu <- as.data.frame(sim.base, out = "qnt", qval = 0.95)

# prevalence
round(data.frame(mean = mn$i.prev[steps],
                 ql = ql$i.prev[steps], qu = qu$i.prev[steps]), 3)

# incidence
round(data.frame(
  mean = (mean(tail(mn$incid, 100)) /
            ((1 - mean(tail(mn$i.prev, 100))) * mean(tail(mn$num, 100)))) * 52 * 100,
  ql = (mean(tail(ql$incid, 100)) /
          ((1 - mean(tail(ql$i.prev, 100))) * mean(tail(ql$num, 100)))) * 52 * 100,
  qu = (mean(tail(qu$incid, 100)) /
          ((1 - mean(tail(qu$i.prev, 100))) * mean(tail(qu$num, 100)))) * 52 * 100), 2)

ir.base <- (sum(mn$incid)/sum((1 - mn$i.prev) * mn$num)) * 52 * 1e5
ir.base

incid.base <- sum(mn$incid)
incid.base


## No Risk Comp
load("data/sim.n2100.rda")
epi_stats(sim, ir.base, incid.base)

## 100% risk comp, all MSM
load("data/sim.n2110.rda")
epi_stats(sim, ir.base, incid.base)

## 100% risk comp, high adherent MSM only
load("data/sim.n2121.rda")
epi_stats(sim, ir.base, incid.base)

boxplot(incid ~ rcomp*hadr, data = df, outline = FALSE)

# Poisson model -----------------------------------------------------------

fn <- list.files("data/", full.names = TRUE)
df <- data.frame(rcomp = NA, hadr = NA, incid = NA, offst = NA)

for (i in 2:length(fn)) {
  load(fn[i])
  rcomp <- rep(sim$param$rcomp.prob, sim$control$nsims)
  hadr <- rep(as.numeric(sim$param$rcomp.hadhr.only), sim$control$nsims)
  sim <- truncate_sim(sim, at = 2600)$epi
  incid <- unname(colSums(sim$incid))
  offst <- unname(colSums((1 - sim$i.prev) * sim$num))
  dft <- data.frame(rcomp, hadr, incid, offst)
  df <- rbind(df, dft)
}

df <- df[-1, ]
df2 <- df[sample(nrow(df), 200), ]

mod <- glm(incid ~ rcomp + rcomp*hadr + offset(log(offst)), family = "poisson", data = df)
summary(mod)

nd0 <- expand.grid(rcomp = seq(0, 1, 0.05), hadr = 0, offst = 5200)
nd1 <- expand.grid(rcomp = seq(0, 1, 0.05), hadr = 1, offst = 5200)
pred0 <- predict(mod, newdata = nd0, type = "response")
pred1 <- predict(mod, newdata = nd1, type = "response")

pred0.pia <- ((ir.base/1000) - pred0)/(ir.base/1000)
pred1.pia <- ((ir.base/1000) - pred1)/(ir.base/1000)

par(mar = c(3,3,2,1), mgp = c(2,1,0))
plot(x = seq(0, 1, 0.05), y = pred0.pia, type = "o", lwd = 2, pch = 20, ylim = c(0, 0.4),
     main = "Percent of Infections Averted by Risk Compensation", col = "blue",
     ylab = "Percent of Infections Averted", xlab = "Relative Increase in Risk Compensation")
lines(x = seq(0, 1, 0.05), y = pred1.pia, col = "red", lwd = 2, type = "o", pch = 20)
abline(h = pred0.pia[1], lty = 2)
legend("bottomleft", legend = c("Baseline PIA", "RC in All", "RC in Highly-Adherent Only"), lty = c(2, 1, 1),
       col = c("black", "blue", "red"), lwd = 2, cex = 0.8, bty = "n")

pred0.pia / pred0.pia[1]
