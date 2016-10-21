
## STI PrEP Figures 1 and 2

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")


# Process Data --------------------------------------------------------

sims <- 1000:1098
cov <- rep(NA, length(sims))
rc <- rep(NA, length(sims))
hr.gc <- rep(NA, length(sims))
hr.ct <- rep(NA, length(sims))
pia.gc <- rep(NA, length(sims))
nnt.gc <- rep(NA, length(sims))
pia.ct <- rep(NA, length(sims))
nnt.ct <- rep(NA, length(sims))

df <- data.frame(sims, cov, rc, hr.gc, hr.ct, pia.gc, pia.ct, nnt.gc, nnt.ct)

load("data/sim.n100.rda")
sim.base <- sim
haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))

for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  df$cov[i] <- sim$param$prep.coverage
  df$rc[i] <- sim$param$rcomp.prob

  # HR
  num.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  denom.gc <- unname(colMeans(tail(sim.base$epi$ir100.gc, 52)))
  vec.hr.gc <- num.gc/denom.gc
  vec.hr.gc <- vec.hr.gc[vec.hr.gc < Inf]
  df$hr.gc[i] <- median(vec.hr.gc, na.rm = TRUE)

  num.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
  denom.ct <- unname(colMeans(tail(sim.base$epi$ir100.ct, 52)))
  vec.hr.ct <- num.ct/denom.ct
  vec.hr.ct <- vec.hr.ct[vec.hr.ct < Inf]
  df$hr.ct[i] <- median(vec.hr.ct, na.rm = TRUE)

  # PIA
  ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
  vec.pia.gc <- vec.nia.gc/ir.base.gc
  vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
  df$pia.gc[i] <- median(vec.pia.gc, na.rm = TRUE)

  ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
  vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
  vec.pia.ct <- vec.nia.ct/ir.base.ct
  vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
  df$pia.ct[i] <- median(vec.pia.ct, na.rm = TRUE)

  # NNT
  py.on.prep <- unname(colSums(sim$epi$prepCurr))/52
  vec.nnt.gc <- py.on.prep / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
  vec.nnt.ct <- py.on.prep / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))

  df$nnt.gc[i] <- median(vec.nnt.gc, na.rm = TRUE)
  df$nnt.ct[i] <- median(vec.nnt.ct, na.rm = TRUE)
  cat("*")
}
df


# Figure 1: PIA by Coverage x Risk Compensation -----------------------

require(gridExtra)
require(lattice)
library(viridis)

pia.loess.gc <- loess(pia.gc ~ cov * rc, data = df, degree = 2, span = 0.15)
pia.fit.gc <- expand.grid(list(cov = seq(0.1, 0.9, 0.002),
                               rc = seq(0, 1, 0.002)))
pia.fit.gc$pia <- as.numeric(predict(pia.loess.gc, newdata = pia.fit.gc))

pia.loess.ct <- loess(pia.ct ~ cov * rc, data = df, degree = 2, span = 0.15)
pia.fit.ct <- expand.grid(list(cov = seq(0.1, 0.9, 0.002),
                               rc = seq(0, 1, 0.002)))
pia.fit.ct$pia <- as.numeric(predict(pia.loess.ct, newdata = pia.fit.ct))
pal <- viridis(n = 21, option = "D")

plot.left <- contourplot(pia ~ cov * rc, data = pia.fit.gc,
                         cuts = 15, region = TRUE,
                         ylab = "Risk Compensation",
                         xlab = "Coverage",
                         main = "Percent GC Infections Averted",
                         col.regions = pal,
                         labels = FALSE,
                         contour = TRUE)

plot.right <- contourplot(pia ~ cov * rc, data = pia.fit.ct,
                          cuts = 15, region = TRUE,
                          ylab = "Risk Compensation",
                          xlab = "Coverage",
                          main = "Percent CT Infections Averted",
                          col.regions = pal,
                          labels = FALSE,
                          contour = TRUE)

tiff(filename = "analysis/Fig1.tiff", height = 6, width = 11, units = "in", res = 250)
grid.arrange(plot.left, plot.right, ncol = 2)
dev.off()


# Figure 2: NNT by Coverage x Risk Compensation -----------------------

nnt.loess.gc <- loess(nnt.gc ~ cov * rc, data = df, degree = 2, span = 0.25)
nnt.fit.gc <- expand.grid(list(cov = seq(0.1, 0.9, 0.002),
                               rc = seq(0, 1, 0.002)))
nnt.fit.gc$nnt <- as.numeric(predict(nnt.loess.gc, newdata = nnt.fit.gc))

nnt.loess.ct <- loess(nnt.ct ~ cov * rc, data = df, degree = 2, span = 0.25)
nnt.fit.ct <- expand.grid(list(cov = seq(0.1, 0.9, 0.002),
                               rc = seq(0, 1, 0.002)))
nnt.fit.ct$nnt <- as.numeric(predict(nnt.loess.ct, newdata = nnt.fit.ct))
pal <- viridis(n = 16, option = "C")

plot.left <- contourplot(nnt ~ cov * rc, data = nnt.fit.gc,
                         cuts = 12, region = TRUE,
                         ylab = "Risk Compensation",
                         xlab = "Coverage",
                         main = "Number Needed to Treat (GC)",
                         col.regions = pal,
                         labels = FALSE)

plot.right <- contourplot(nnt ~ cov * rc, data = nnt.fit.ct,
                          cuts = 12, region = TRUE,
                          ylab = "Risk Compensation",
                          xlab = "Coverage",
                          main = "Number Needed to Treat (CT)",
                          col.regions = pal,
                          labels = FALSE)

tiff(filename = "analysis/Fig2.tiff", height = 6, width = 11, units = "in", res = 250)
grid.arrange(plot.left, plot.right, ncol = 2)
dev.off()

library(plotly)
plot_ly(x = pia.fit.gc$cov, y = pia.fit.gc$rc, z = pia.fit.gc$pia, type = "contour")


# Supp Figure: HR by Coverage x Risk Compensation ---------------------

hr.loess.gc <- loess(hr.gc ~ cov * rc, data = df, degree = 2, span = 0.2)
hr.fit.gc <- expand.grid(list(cov = seq(0.1, 0.9, 0.002),
                              rc = seq(0, 1, 0.002)))
hr.fit.gc$hr <- as.numeric(predict(hr.loess.gc, newdata = hr.fit.gc))

hr.loess.ct <- loess(hr.ct ~ cov * rc, data = df, degree = 2, span = 0.2)
hr.fit.ct <- expand.grid(list(cov = seq(0.1, 0.9, 0.002),
                              rc = seq(0, 1, 0.002)))
hr.fit.ct$hr <- as.numeric(predict(hr.loess.ct, newdata = hr.fit.ct))
pal <- viridis(n = 12, option = "C")

plot.left <- contourplot(hr ~ cov * rc, data = hr.fit.gc,
                         cuts = 9, region = TRUE,
                         ylab = "Risk Compensation",
                         xlab = "Coverage",
                         main = "GC Hazard Ratio",
                         col.regions = pal,
                         labels = FALSE)

plot.right <- contourplot(hr ~ cov * rc, data = hr.fit.ct,
                          cuts = 9, region = TRUE,
                          ylab = "Risk Compensation",
                          xlab = "Coverage",
                          main = "CT Hazard Ratio",
                          col.regions = pal,
                          labels = FALSE)

tiff(filename = "analysis/Fig1b.tiff", height = 6, width = 11, units = "in", res = 250)
grid.arrange(plot.left, plot.right, ncol = 2)
dev.off()
