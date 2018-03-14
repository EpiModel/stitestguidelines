rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
#source("analysis/fx.R")


# 4000 - Base
# 4001 - Base + 5% ann
# 4002 - Base + 10% ann
# 4003 - Base + 15% ann
# 4004 - Base + 20% ann
# 4005 - Base + 25% ann
# 4006 - Base + 30% ann
# 4007 - Base + 35% ann
# 4008 - Base + 40% ann

# 4009 - Base + 5% HR
# 4018 - Base + 10% HR
# 4027 - Base + 15% HR
# 4036 - Base + 20% HR
# 4045 - Base + 25% HR
# 4054 - Base + 30% HR
# 4063 - Base + 35% HR
# 4072 - Base + 40% HR
# 4081 - Base + 45% HR
# 4090 - Base + 50% HR
# 4099 - Base + 55% HR
# 4108 - Base + 60% HR
# 4117 - Base + 65% HR
# 4126 - Base + 70% HR
# 4135 - Base + 75% HR
# 4144 - Base + 80% HR
# 4153 - Base + 85% HR
# 4162 - Base + 90% HR
# 4171 - Base + 95% HR
# 4180 - Base + 100% HR

## Base: n4000
## Varying Ann cov: n4001 - 4008
## Varying HR cov: 4009, 4018, 4027, 4036, 4045, 4054, 4063, 4072, 4081, 4090, 4099, 4108, 4117, 4126, 4135, 4144, 4153, 4162, 4171, 4180

tiff(filename = "data/AnnIR.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- 3000:3008
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "NG/CT Incidence by Annual Screening Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("Base", "+5%", "+10%", "+15%", "+20%",
                                "+25%", "+30%", "+35%", "+40%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

sims <- 3000:3008
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.syph", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "Syph Incidence by Annual Screening Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("Base", "+5%", "+10%", "+15%", "+20%",
                                "+25%", "+30%", "+35%", "+40%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()


tiff(filename = "data/HRIR.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(4000, 4009, 4018, 4027, 4036, 4045, 4054, 4063, 4072, 4081, 4090, 4099,
          4108, 4117, 4126, 4135, 4144)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "NG/CT Incidence by HR Screening Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("Base", "+5%", "+10%", "+15%", "+20%",
                                "+25%", "+30%", "+35%", "+40%", "+45%",
                                "+50%", "+55%", "+60%", "+65%", "+70%",
                                "+75%", "+80%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

sims <- c(4000, 4009, 4018, 4027, 4036, 4045, 4054, 4063, 4072, 4081, 4090, 4099,
          4108, 4117, 4126, 4135, 4144)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.syph", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
       main = "Syph Incidence by HR Screening Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 6))
}
legend("bottomleft", legend = c("Base", "+5%", "+10%", "+15%", "+20%",
                                "+25%", "+30%", "+35%", "+40%", "+45%",
                                "+50%", "+55%", "+60%", "+65%", "+70%",
                                "+75%", "+80%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()
