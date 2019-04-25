# Process Data

## Supplemental Figure 2 --------------------------------------------------------
## U-Shaped Diagram
## Line plot: Number of partners and PIA
## Partner cutoff by HR coverage
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")
library("wesanderson")

# Baseline
load("data/followup/Guidelines Paper/sim.n9000.rda")
sim.base <- sim

# Box Plots of partner number and PIA
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
mn.base <- as.data.frame(sim.base)
# Note, this doesn't account for dually-infected
ir.base <- ((colSums(sim.base$epi$incid.rgc, na.rm = TRUE) +
              colSums(sim.base$epi$incid.ugc, na.rm = TRUE) +
              colSums(sim.base$epi$incid.rct, na.rm = TRUE) +
              colSums(sim.base$epi$incid.uct, na.rm = TRUE)) /
              sum((1 - mn.base$prev.gc - mn.base$prev.ct) * mn.base$num)) * 52 * 1e5
incid.base.sti <- (colSums(sim.base$epi$incid.rgc, na.rm = TRUE) +
                     colSums(sim.base$epi$incid.ugc, na.rm = TRUE) +
                     colSums(sim.base$epi$incid.rct, na.rm = TRUE) +
                     colSums(sim.base$epi$incid.uct, na.rm = TRUE))
tests.base.sti <- (colSums(sim.base$epi$GCasympttests, na.rm = TRUE) +
                     colSums(sim.base$epi$CTasympttests, na.rm = TRUE))

# Partner Number threshold:
sims <- c(9000, 9009, 9037:9040) #, 9037:9045)
df.pia <- data.frame(rep(NA, 256))
df.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {

  load(list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE))
  mn <- as.data.frame(sim)
  ir <- (colSums(sim$epi$incid.rgc, na.rm = TRUE) +
           colSums(sim$epi$incid.ugc, na.rm = TRUE) +
           colSums(sim$epi$incid.rct, na.rm = TRUE) +
           colSums(sim$epi$incid.uct, na.rm = TRUE)) /
    sum((1 - mn$prev.gc - mn$prev.ct)  * mn$num) * 52 * 1e5 #Note that this doesn't account for dually-infected
  vec.nia <- round(ir.base - unname(ir), 1)
  df.pia[, i] <- vec.nia / ir.base

  # stiasympttests <- sum(mn$GCasympttests) + sum(mn$CTasympttests)
  # incid.sti <- sum(mn$incid.rgc) + sum(mn$incid.ugc) + sum(mn$incid.rct) + sum(mn$incid.uct)
  incid.sti <- (colSums(sim$epi$incid.rgc, na.rm = TRUE) +
                  colSums(sim$epi$incid.ugc, na.rm = TRUE) +
                  colSums(sim$epi$incid.rct, na.rm = TRUE) +
                  colSums(sim$epi$incid.uct, na.rm = TRUE))
  stiasympttests <- (colSums(sim$epi$GCasympttests, na.rm = TRUE) +
                       colSums(sim$epi$CTasympttests, na.rm = TRUE))
  df.nnt[, i] <- (stiasympttests - tests.base.sti) / (incid.base.sti - incid.sti)

}

names(df.pia) <- names(df.nnt) <- c("No HR Screening", ">1",
                                    ">2", ">3",
                                    ">4", ">5")#, ">6 partners",
                                    # ">7 partners", ">8 partners", ">9 partners",
                                    # ">10 partners")
head(df.pia)
head(df.nnt)

# Set comparator column to zero
df.nnt[is.nan(df.nnt$`No HR Screening`), 1] <- 0

pal <- wes_palette("Zissou1")[c(1, 5)]
tiff(filename = "analysis/Supp Fig 2.tiff", height = 4, width = 8, units = "in", res = 250)

par(mfrow = c(1, 2), mar = c(3,3,2.5,1), oma = c(0, 0, 3, 0), mgp = c(2,1,0))

# Left Panel: PIA
boxplot(df.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "PIA by High-Risk Partner Threshold",
        xlab = "Partner Threshold", ylab = "Percent of Infections Averted")

# Right Panel: NNT
boxplot(df.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "NNT by High-Risk Partner Threshold",
        xlab = "Partner Threshold", ylab = "Number Needed to Screen")
mtext("Supplementary Figure 2: Boxplots of Percent of Infections Averted (PIA)
      and Number Needed to Screen (NNS)", outer = TRUE, cex = 1)
dev.off()

## Supplemental Figure 3 -------------------------------------------------------
## Line plot: Incidence and SA interval/HR Interval
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# 9009 - 5% HR, 1 partner
## SA Screening Interval: 9029:9032
tiff(filename = "analysis/Supp Fig 3.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(3,3,2,1.2), oma = c(0, 0, 2, 0),mgp = c(2,1,0))
sims <- c(9029:9030, 9009, 9031:9032)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.1,
       main = "STI Incidence by Lower-Risk Screening Interval",
       xlab = "Week", ylab = "Incidence Rate (IR) per 100 Person-Years at risk", ylim = c(0, 13))
}
legend("bottomleft", legend = c("6 months", "9 months", "12 months",
                                "15 months", "18 months"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

## HR Screening Interval: 9033:9036
## 5% HR, 1 partner, 9009
sims <- c(9033:9034, 9009, 9035:9036)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.1,
       main = "STI Incidence by Higher-Risk Screening Interval",
       xlab = "Week", ylab = "Incidence Rate (IR) per 100 Person-Years at risk", ylim = c(0, 13))
}
legend("bottomleft", legend = c("1 months", "3 months", "6 months",
                                "9 months", "12 months"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
mtext("Supplementary Figure 3: STI Incidence by STI Screening Interval",
      outer = TRUE, cex = 1.5)
dev.off()

## Supplemental Figure 4 -------------------------------------------------------
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
tiff(filename = "analysis/Supp Fig 4.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), oma = c(0, 0, 2, 0), mgp = c(2,1,0))
sims <- c(9000, 9002, 9004, 9006, 9008)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.1,
       main = "STI Incidence by Lower-Risk Screening Coverage",
       xlab = "Week", ylab = "Incidence Rate (IR) per 100 Person-Years at risk", ylim = c(0, 13))
}
legend("bottomleft", legend = c("Baseline", "10% increase", "20% increase",
                                "30% increase", "40% increase"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")


## Varying HR coverage: 9009:9028
## 20% = 9012, 40% = 9016, 60% = 9020, 80% = 9024, 100% = 9028
sims <- c(9000, 9012, 9016, 9020, 9024, 9028)
# sims <- c(9000, 9009, 9010, 9012, 9014, 9016, 9018, 9020, 9022, 9024, 9026, 9028)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.1,
       main = "STI Incidence by High-Risk Screening Coverage",
       xlab = "Week", ylab = "Incidence Rate (IR) per 100 Person-Years at risk", ylim = c(0, 13))
}
legend(x = 330, y = 9, legend = c("0%", "20%", "40%", "60%", "80%", "100%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
mtext("Supplementary Figure 4: STI Incidence by STI Screening Coverage",
      outer = TRUE, cex = 1.5)
dev.off()


## Supplemental Figure 5 -------------------------------------------------------
## Line plot: Incidence and partner cutoff
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# Baseline for this (1 partner): 9000
## Varying Partner Cutoffs: 9037:9045
tiff(filename = "analysis/Supp Fig 5.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,1), mar = c(3,3,2,1.2), oma = c(0, 0, 2, 0), mgp = c(2,1,0))
sims <- c(9000, 9009, 9037:9040)#9037:9045)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.1,
       #main = "Supplementary Figure 5: STI Incidence by High-Risk Partner Threshold",
       xlab = "Week", ylab = "Incidence Rate (IR) per 100 Person-Years at risk", ylim = c(0, 13))
}
legend("bottomleft", legend = c("No HR Screening",
                                ">1 partner", ">2 partners", ">3 partners",
                                ">4 partners", ">5 partners"), #, ">6 partners",
                                # ">7 partners", ">8 partners", ">9 partners",
                                # ">10 partners"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
mtext("Supplementary Figure 5: STI Incidence by High-Risk Partner Threshold",
      outer = TRUE, cex = 1.5)
dev.off()


## Supplemental Figure 6 -------------------------------------------------------
## NNS
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")


## Statistic on proportion of screening tests ----------------------------------
load("data/followup/Guidelines Paper/sim.n9000.rda")

asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE)) + unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
sympt.tests <- unname(colSums(sim$epi$CTsympttests, na.rm = TRUE)) + unname(colSums(sim$epi$GCsympttests, na.rm = TRUE))
frac <- asympt.tests / (asympt.tests + sympt.tests)

quantile(frac, probs = c(0.025, 0.5, 0.975)) #94.38 (91.87 - 96.46)
