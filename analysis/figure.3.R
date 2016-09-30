

## STI PrEP Figure 3

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Varying STI testing interval: n2000 to 2004
## Varying prob STI treatment: n2005 to 2009

tiff(filename = "analysis/Fig3.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- 2000:2004
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.sti", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "Combined GC/CT Incidence by STI Testing Interval",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("1 mo", "3 mo", "6 mo", "9 mo", "12 mo"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

sims <- 2005:2009
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.sti", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "Combined GC/CT Incidence by Proportion Screened",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 15))
}
legend("bottomleft", legend = c("0%", "25%", "50%", "75%", "100%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()
