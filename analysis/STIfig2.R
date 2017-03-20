## STI Testing guidelines Figures 1 and 2

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")


# Process Data --------------------------------------------------------

sims <- c(3011:3130)
anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
hr.gc <- rep(NA, length(sims))
hr.ct <- rep(NA, length(sims))
hr.syph <- rep(NA, length(sims))
pia.gc <- rep(NA, length(sims))
nnt.gc <- rep(NA, length(sims))
pia.ct <- rep(NA, length(sims))
nnt.ct <- rep(NA, length(sims))
pia.syph <- rep(NA, length(sims))
nnt.syph <- rep(NA, length(sims))

df <- data.frame(sims, anncov, hrcov, hr.gc, hr.ct, hr.syph, pia.gc, pia.ct, nnt.gc, nnt.ct)

load("data/followup/sim.n3000.rda")
sim.base <- sim
haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    df$cov[i] <- sim$param$stianntest.coverage
    df$rc[i] <- sim$param$stihighrisktest.coverage
    
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
    
    num.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
    denom.syph <- unname(colMeans(tail(sim.base$epi$ir100.syph, 52)))
    vec.hr.syph <- num.syph/denom.syph
    vec.hr.syph <- vec.hr.syph[vec.hr.syph < Inf]
    df$hr.syph[i] <- median(vec.hr.syph, na.rm = TRUE)
    
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
    
    ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
    vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
    vec.pia.syph <- vec.nia.syph/ir.base.syph
    vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
    df$pia.syph[i] <- median(vec.pia.syph, na.rm = TRUE)
    
    # NNT
    gc.tests <- unname(tail(sim$epi$totalGCasympttests, 1))
    ct.tests <- unname(tail(sim$epi$totalCTasympttests, 1))
    syph.tests <- unname(tail(sim$epi$totalsyphasympttests, 1))
    vec.nnt.gc <- gc.tests / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
    vec.nnt.ct <- ct.tests / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))
    vec.nnt.syph <- syph.tests / (median(incid.base.syph) - unname(colSums(sim$epi$incid.syph)))
    
    df$nnt.gc[i] <- median(vec.nnt.gc, na.rm = TRUE)
    df$nnt.ct[i] <- median(vec.nnt.ct, na.rm = TRUE)
    df$nnt.syph[i] <- median(vec.nnt.syph, na.rm = TRUE)
    cat("*")
}
df


# Figure 2a: PIA by Coverage x Risk Compensation -----------------------

require(gridExtra)
require(lattice)
library(viridis)

pia.loess.gc <- loess(pia.gc ~ hrcov * anncov, data = df, degree = 2, span = 0.15)
pia.fit.gc <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                               hrcov = seq(0, 1, 0.002)))
pia.fit.gc$pia <- as.numeric(predict(pia.loess.gc, newdata = pia.fit.gc))

pia.loess.ct <- loess(pia.ct ~ hrcov * anncov, data = df, degree = 2, span = 0.15)
pia.fit.ct <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                                hrcov = seq(0, 1, 0.002)))
pia.fit.ct$pia <- as.numeric(predict(pia.loess.ct, newdata = pia.fit.ct))

pia.loess.syph <- loess(pia.syph ~ hrcov * anncov, data = df, degree = 2, span = 0.15)
pia.fit.syph <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                               hrcov = seq(0, 1, 0.002)))
pia.fit.syph$pia <- as.numeric(predict(pia.loess.syph, newdata = pia.fit.syph))

pal <- viridis(n = 21, option = "D")

plot.left <- contourplot(pia ~ hrcov * anncov, data = pia.fit.gc,
                         cuts = 15, region = TRUE,
                         ylab = "Higher-Risk Testing Coverage",
                         xlab = "Lower-Risk Testing Coverage",
                         main = "Percent GC Infections Averted",
                         col.regions = pal,
                         labels = FALSE,
                         contour = TRUE)

plot.center <- contourplot(pia ~ hrcov * anncov, data = pia.fit.ct,
                          cuts = 15, region = TRUE,
                          ylab = "Higher-Risk Testing Coverage",
                          xlab = "Lower-Risk Testing Coverage",
                          main = "Percent CT Infections Averted",
                          col.regions = pal,
                          labels = FALSE,
                          contour = TRUE)

plot.right <- contourplot(pia ~ hrcov * anncov, data = pia.fit.syph,
                           cuts = 15, region = TRUE,
                           ylab = "Higher-Risk Testing Coverage",
                           xlab = "Lower-Risk Testing Coverage",
                           main = "Percent Syph Infections Averted",
                           col.regions = pal,
                           labels = FALSE,
                           contour = TRUE)

tiff(filename = "analysis/Fig2a.tiff", height = 6, width = 11, units = "in", res = 250)
grid.arrange(plot.left, plot.right, ncol = 3)
dev.off()


# Figure 2b: NNT by Coverage x Risk Compensation -----------------------

nnt.loess.gc <- loess(nnt.gc ~ hrcov * anncov, data = df, degree = 2, span = 0.25)
nnt.fit.gc <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                               hrcov = seq(0, 1, 0.002)))
nnt.fit.gc$nnt <- as.numeric(predict(nnt.loess.gc, newdata = nnt.fit.gc))

nnt.loess.ct <- loess(nnt.ct ~ hrcov * anncov, data = df, degree = 2, span = 0.25)
nnt.fit.ct <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                               hrcov = seq(0, 1, 0.002)))
nnt.fit.ct$nnt <- as.numeric(predict(nnt.loess.ct, newdata = nnt.fit.ct))

nnt.loess.syph <- loess(nnt.syph ~ hrcov * anncov, data = df, degree = 2, span = 0.25)
nnt.fit.syph <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                               hrcov = seq(0, 1, 0.002)))
nnt.fit.syph$nnt <- as.numeric(predict(nnt.loess.syph, newdata = nnt.fit.syph))

pal <- viridis(n = 16, option = "C")

plot.left <- contourplot(nnt ~ hrcov * anncov, data = nnt.fit.gc,
                         cuts = 12, region = TRUE,
                         ylab = "Higher-Risk Testing Coverage",
                         xlab = "Lower-Risk Testing Coverage",
                         main = "Number Needed to Treat (GC)",
                         col.regions = pal,
                         labels = FALSE)

plot.center <- contourplot(nnt ~ hrcov * anncov, data = nnt.fit.ct,
                          cuts = 12, region = TRUE,
                          ylab = "Higher-Risk Testing Coverage",
                          xlab = "Lower-Risk Testing Coverage",
                          main = "Number Needed to Treat (CT)",
                          col.regions = pal,
                          labels = FALSE)

plot.center <- contourplot(nnt ~ hrcov * anncov, data = nnt.fit.syph,
                           cuts = 12, region = TRUE,
                           ylab = "Higher-Risk Testing Coverage",
                           xlab = "Lower-Risk Testing Coverage",
                           main = "Number Needed to Treat (syph)",
                           col.regions = pal,
                           labels = FALSE)

tiff(filename = "stiguidelines/analysis/Fig2b.tiff", height = 6, width = 11, units = "in", res = 250)
grid.arrange(plot.left, plot.right, ncol = 3)
dev.off()

library(plotly)
plot_ly(x = pia.fit.gc$anncov, y = pia.fit.gc$hrcov, z = pia.fit.gc$pia, type = "contour")
plot_ly(x = pia.fit.ct$anncov, y = pia.fit.ct$hrcov, z = pia.fit.ct$pia, type = "contour")
plot_ly(x = pia.fit.syph$anncov, y = pia.fit.syph$hrcov, z = pia.fit.syph$pia, type = "contour")



# Supp Figure: HR by Coverage x Risk Compensation ---------------------

hr.loess.gc <- loess(hr.gc ~ hrcov * anncov, data = df, degree = 2, span = 0.2)
hr.fit.gc <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                              hrcov = seq(0, 1, 0.002)))
hr.fit.gc$hr <- as.numeric(predict(hr.loess.gc, newdata = hr.fit.gc))

hr.loess.ct <- loess(hr.ct ~ hrcov * anncov, data = df, degree = 2, span = 0.2)
hr.fit.ct <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                              hrcov = seq(0, 1, 0.002)))
hr.fit.ct$hr <- as.numeric(predict(hr.loess.ct, newdata = hr.fit.ct))

hr.loess.syph <- loess(hr.syph ~ hrcov * anncov, data = df, degree = 2, span = 0.2)
hr.fit.syph <- expand.grid(list(anncov = seq(0.0, 1.0, 0.002),
                              hrcov = seq(0, 1, 0.002)))
hr.fit.syph$hr <- as.numeric(predict(hr.loess.syph, newdata = hr.fit.syph))

pal <- viridis(n = 12, option = "C")

plot.left <- contourplot(hr ~ hrcov * anncov, data = hr.fit.gc,
                         cuts = 9, region = TRUE,
                         ylab = "Higher-Risk Testing Coverage",
                         xlab = "Lower-Risk Testing Coverage",
                         main = "GC Hazard Ratio",
                         col.regions = pal,
                         labels = FALSE)

plot.center <- contourplot(hr ~ hrcov * anncov, data = hr.fit.ct,
                          cuts = 9, region = TRUE,
                          ylab = "Higher-Risk Testing Coverage",
                          xlab = "Lower-Risk Testing Coverage",
                          main = "CT Hazard Ratio",
                          col.regions = pal,
                          labels = FALSE)

plot.right <- contourplot(hr ~ hrcov * anncov, data = hr.fit.syph,
                          cuts = 9, region = TRUE,
                          ylab = "Higher-Risk Testing Coverage",
                          xlab = "Lower-Risk Testing Coverage",
                          main = "Syph Hazard Ratio",
                          col.regions = pal,
                          labels = FALSE)


tiff(filename = "analysis/Fig2b.tiff", height = 6, width = 11, units = "in", res = 250)
grid.arrange(plot.left, plot.right, ncol = 3)
dev.off()