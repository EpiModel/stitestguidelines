# Process Data

## Supplemental Figure 1 --------------------------------------------------------
## U-Shaped Diagram
## Line plot: Number of partners and PIA
## Partner cutoff by HR coverage
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# Baseline
load("data/followup/Guidelines Paper/sim.n9000.rda")
sim.base <- sim

# Box Plots of partner number and PIA
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
mn.base <- as.data.frame(sim.base)
ir.base <- (sum(mn.base$incid.sti)/sum((1 - mn.base$prev.sti) * mn.base$num)) * 52 * 1e5
incid.base.sti <- sum(mn.base$incid.sti)
tests.base.sti <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))

# incid.sti <- unname(colSums(sim$epi$incid.sti, na.rm = TRUE))
# vec.nia.ct <- incid.base.ct - incid.ct
# vec.pia.ct <- vec.nia.ct/incid.base.ct
# pia.ct <- median(vec.pia.ct, na.rm = TRUE)

# Partner Number threshold:
sims <- c(9000, 9037:9045)
df.pia <- data.frame(rep(NA, 256))
df.nnt <- data.frame(rep(NA, 256))
for (i in seq_along(sims)) {
  load(list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE))
  mn <- as.data.frame(sim)
  ir <- (colSums(sim$epi$incid.sti, na.rm = TRUE)) /
    sum((1 - mn$prev.sti)  * mn$num) * 52 * 1e5
  vec.nia <- round(ir.base - unname(ir), 1)
  df.pia[, i] <- vec.nia / ir.base


  stiasympttests <- unname(colSums(sim$epi$stiasympttests))
  df.nnt[, i] <- (stiasympttests - tests.base.sti) / (incid.base.sti - unname(colSums(sim$epi$incid.sti)))

}
names(df.pia) <- names(df.nnt) <- c(">1 partner", ">2 partners", ">3 partners",
                                    ">4 partner", ">5 partners", ">6 partners",
                                    ">7 partner", ">8 partners", ">9 partners",
                                    ">10 partners")
head(df.pia)
head(df.nnt)

pal <- wes_palette("Zissou")[c(1, 5)]
tiff(filename = "analysis/Supp Figure 1", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))

# Left Panel: PIA
boxplot(df.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "PIA by High-Risk Partner Threshold",
        xlab = "Partner Threshold", ylab = "Percent Infections Averted")

# Right Panel: NNT
boxplot(df.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "NNS by High-Risk Partner Threshold",
        xlab = "Partner Threshold", ylab = "Percent Infections Averted")

dev.off()

## Supplemental Figure 2 -------------------------------------------------------
## Line plot: Incidence and SA interval/HR Interval
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# 9009 - 5% HR, 1 partner
## SA Screening Interval: 9029:9032
tiff(filename = "analysis/Supp Fig 2", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(9029:9030, 9009, 9031:9032)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "STI Incidence by Lower-Risk Screening Interval",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("6 mo", "9 mo", "12 mo (Baseline)", "15 mo", "18 mo"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

## HR Screening Interval: 9033:9036
## 5% HR, 1 partner, 9009
sims <- c(9033:9034, 9009, 9035:9036)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "STI Incidence by Higher-Risk Screening Interval",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("1 mo", "3 mo", "6 mo (Baseline)", "9 mo", "12 mo"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()

## Supplemental Figure 3 -------------------------------------------------------
## Line plot: Incidence and coverage
## 2 panels (One for HR and one for SA)
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# Baseline for this (1 partner): 9000
## Varying SA coverage: 9001:9008
tiff(filename = "analysis/Supp Fig 3", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))

sims <- c(9000, 9002, 9004, 9006, 9008)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "STI Incidence by Lower-Risk Screening Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("Baseline", "10% increase", "20% increase",
                                "30% increase", "40% increase"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")


## Varying HR coverage: 9009:9028
sims <- c(9000, 9010, 9012, 9014, 9016, 9018, 9020, 9022, 9024, 9026, 9028)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "STI Incidence by High-Risk Partner Threshold",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("5% (Baseline)", "10%", "20%", "30%", "40%",
                                "50%", "60%", "70%", "80%", "90%", "100%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

dev.off()


## Supplemental Figure 4 -------------------------------------------------------
## Line plot: Incidence and partner cutoff
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# Baseline for this (1 partner): 9000
## Varying Partner Cutoffs: 9037:9045
tiff(filename = "analysis/Supp Fig 4", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(9000, 9037:9045)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "STI Incidence by High-Risk Partner Threshold",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c(">1 partner", ">2 partners", ">3 partners",
                                ">4 partner", ">5 partners", ">6 partners",
                                ">7 partner", ">8 partners", ">9 partners",
                                ">10 partners"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

dev.off()


## Supplemental Figure 5 -------------------------------------------------------
## NNS
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")


## Statistic on proportion of screening tests ----------------------------------
load("data/followup/Guidelines Paper/sim.n9000.rda")

asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))
sympt.tests <- unname(colSums(sim$epi$stisympttests, na.rm = TRUE))
frac <- asympt.tests / (asympt.tests + sympt.tests)

quantile(frac, probs = c(0.025, 0.5, 0.975))
