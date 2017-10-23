# NEEMA
### Testing in the last 12 months ----------------------------------------------

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
load("data/followup/sim.n3000.rda")
sim.base <- sim

# Get denominators for each eligibility indication
quantile(colMeans(tail(sim.base$epi$recentpartners.prop)))
quantile(colMeans(tail(sim.base$epi$stiactiveind.prop)))
quantile(colMeans(sim.base$epi$recentpartners.prop[2:53, ]))
quantile(colMeans(sim.base$epi$stiactiveind.prop[2:53, ]))

# Varying Lower-Risk Coverage
# 3002, 3004, 3006, 3008
# : Annual = 10%, 20%, 30%, 40% increase, 364 days, HR = 0%, 182 days
# Varying Higher-Risk Coverage
# 3018, 3036, 3054, 3072
#  Higher-risk = 0.1 - 0.4, by 0.1, 182 days, Ann = 10%, 364 days

# Screening Intervals:
# 3191 as reference
# 3189, 3190, 3191 (ref), 3192, 3193: Annual = 182 days, 273 days, 364 days (ref), 448 days, 539 days, HR = 50%, ANN = 50%, 182 days

# Newer way:

sims <- c(3000, # Baseline
          3002, 3004, 3006, 3008, # Varying LR coverage
          3018, 3036, 3054, 3072, # Varying HR coverage
          3189, 3190, 3192, 3193 # Varying interval
)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
tst.gc.12mo.hivpos <- rep(NA, length(sims))
tst.ct.12mo.hivpos <- rep(NA, length(sims))
tst.syph.12mo.hivpos <- rep(NA, length(sims))
tst.gc.12mo.hivneg <- rep(NA, length(sims))
tst.ct.12mo.hivneg <- rep(NA, length(sims))
tst.syph.12mo.hivneg <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint,
                 tst.gc.12mo.hivpos, tst.ct.12mo.hivpos, tst.syph.12mo.hivpos,
                 tst.gc.12mo.hivneg, tst.ct.12mo.hivneg, tst.syph.12mo.hivneg)

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int

  # Tested in the last 12 months
  df$tst.gc.12mo.hivneg[i] <- median(unname(colMeans(tail(sim$epi$test.gc.12mo.hivneg, 52))))
  df$tst.ct.12mo.hivneg[i] <- median(unname(colMeans(tail(sim$epi$test.ct.12mo.hivneg, 52))))
  df$tst.syph.12mo.hivneg[i] <- median(unname(colMeans(tail(sim$epi$test.syph.12mo.hivneg, 52))))
  df$tst.gc.12mo.hivpos[i] <- median(unname(colMeans(tail(sim$epi$test.gc.12mo.hivpos, 52))))
  df$tst.ct.12mo.hivpos[i] <- median(unname(colMeans(tail(sim$epi$test.ct.12mo.hivpos, 52))))
  df$tst.syph.12mo.hivpos[i] <- median(unname(colMeans(tail(sim$epi$test.syph.12mo.hivpos, 52))))

  cat("*")

}

df
#write.csv(df, "analysis/NEEMA Testing Table.csv")


## Incidence by testing group --------------------------------------------------
rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Base STI lower-risk testing interval (364 days):
## Varying STI lower-risk testing interval
#tiff(filename = "analysis/Fig3a.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3000) # Choose at 5% HR
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
       lty = 1,
       main = "STI Incidence by Lower-Risk STI Screening Interval",
       xlab = "Week", ylab = "IR per 100 PYAR")
  #abline(h = seq(0, 8, 0.2), lty = 2, col = "gray")
  plot(sim, y = "ir100.sti.tttraj1", add = TRUE, ylim = c(0, 8),
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
       lty = 2,
       main = "STI Incidence by Lower-Risk STI Screening Interval \n 50% Ann / 0% HR Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR")
  plot(sim, y = "ir100.sti.tttraj2", add = TRU, ylim = c(0, 8),
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
       lty = 3,
       main = "STI Incidence by Lower-Risk STI Screening Interval \n 50% Ann / 0% HR Coverage",
       xlab = "Week", ylab = "IR per 100 PYAR")
}
legend("bottomleft", legend = c("182 days", "273 days", "364 days", "448 days", "539 days", "Base - 10% Ann, 364 days"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

#
# ## Base STI higher-risk testing interval
# ## Varying STI higher-risk testing interval
# sims <- c(3447:3451, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
#        main = "STI Incidence by Higher-Risk STI Screening Interval \n 50% HR / 10% Ann Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 8, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("28 days", "91 days", "182 days", "273 days", "364 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()
#
#
# ### HIV
# ## Base STI lower-risk testing interval (364 days):
# ## Varying STI lower-risk testing interval:
# tiff(filename = "analysis/Fig3b.tiff", height = 6, width = 11, units = "in", res = 250)
# par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
# sims <- c(3442:3446, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
#
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "HIV Incidence by Lower-Risk STI Screening Interval \n 50% Ann / 0% HR Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 4, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("182 days", "273 days", "364 days", "448 days", "539 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
#
#
# ## Base STI higher-risk testing interval: 3003 - 10% annual screening
# ## Varying STI higher-risk testing interval: 3447:3451
# sims <- c(3447:3451, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "HIV Incidence by Higher-Risk STI Screening Interval \n 50% HR / 10% Ann Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 4, 0.2), lty = 2, col = "gray")
# }
#
# legend("bottomleft", legend = c("28 days", "91 days", "182 days", "273 days", "364 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()













# ## STI Test Guidelines Figures
#
# rm(list = ls())
# suppressMessages(library("EpiModelHIV"))
# library("EpiModelHPC")
# library("dplyr")
# source("analysis/fx.R")
#
# ## Treatment Progression -------------------------------------------------------
#
# rm(list = ls())
# library("EpiModelHIV")
# library("EpiModelHPC")
# library("dplyr")
# source("analysis/fx.R")
#
# ## Base STI treatment completion: 3003
# ## Varying Treatment completion for 50% Annual
# tiff(filename = "analysis/Fig4a.tiff", height = 6, width = 11, units = "in", res = 250)
# par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
# sims <- c(3452:3472, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
#
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "STI Incidence by Treatment Probability\n 50% Ann / 0% HR Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(4, 8, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("0%","5%", "10%", "15%",
#                                 "20%", "25%", "30%", "35%",
#                                 "40%", "45%", "50%", "55%",
#                                 "60%", "65%", "70%", "75%",
#                                 "80%", "85%", "90%", "95%", "100%", "Base - 10% Ann"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
#
# ## Varying Treatment completion for 10% higher-risk and 50% lower-risk
# sims <- c(3473:3493, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "STI Incidence by Treatment Probability \n 50% HR / 10% Ann Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(4, 8, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("0%","5%", "10%", "15%",
#                                 "20%", "25%", "30%", "35%",
#                                 "40%", "45%", "50%", "55%",
#                                 "60%", "65%", "70%", "75%",
#                                 "80%", "85%", "90%", "95%", "Base - 10% Ann"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()
#
#
# ##### HIV
#
# ## Varying Treatment completion for 50% lower risk and 0% higher-risk: n3175-n3184
# tiff(filename = "analysis/Fig4b.tiff", height = 6, width = 11, units = "in", res = 250)
# par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
# sims <- c(3452:3472, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
#
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "HIV Incidence by STI Treatment Probability \n 50% Ann / 0% HR Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 4, 0.1), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("0%","5%", "10%", "15%",
#                                 "20%", "25%", "30%", "35%",
#                                 "40%", "45%", "50%", "55%",
#                                 "60%", "65%", "70%", "75%",
#                                 "80%", "85%", "90%", "95%", "Base - 10% Ann"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
#
# ## Varying Treatment completion for 50% higher-risk and 10% lower-risk: n3185-n3194
# sims <- c(3473:3493, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "HIV Incidence by STI Treatment Probability \n 50% HR / 10% Ann Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 4, 0.1), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("0%","5%", "10%", "15%",
#                                 "20%", "25%", "30%", "35%",
#                                 "40%", "45%", "50%", "55%",
#                                 "60%", "65%", "70%", "75%",
#                                 "80%", "85%", "90%", "95%", "Base - 10% Ann"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()
#
#
#
# ## Screening Intervals -------------------------------------------------------
#
# rm(list = ls())
# library("EpiModelHIV")
# library("EpiModelHPC")
# library("dplyr")
# source("analysis/fx.R")
#
# ## Base STI lower-risk testing interval (364 days):
# ## Varying STI lower-risk testing interval
# tiff(filename = "analysis/Fig3a.tiff", height = 6, width = 11, units = "in", res = 250)
# par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
# sims <- c(3442:3446, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
#
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
#        main = "STI Incidence by Lower-Risk STI Screening Interval \n 50% Ann / 0% HR Coverage",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 8, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("182 days", "273 days", "364 days", "448 days", "539 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
#
#
# ## Base STI higher-risk testing interval
# ## Varying STI higher-risk testing interval
# sims <- c(3447:3451, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
#        main = "STI Incidence by Higher-Risk STI Screening Interval \n 50% HR / 10% Ann Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 8, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("28 days", "91 days", "182 days", "273 days", "364 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()
#
#
# ### HIV
# ## Base STI lower-risk testing interval (364 days):
# ## Varying STI lower-risk testing interval:
# tiff(filename = "analysis/Fig3b.tiff", height = 6, width = 11, units = "in", res = 250)
# par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
# sims <- c(3442:3446, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
#
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "HIV Incidence by Lower-Risk STI Screening Interval \n 50% Ann / 0% HR Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 4, 0.2), lty = 2, col = "gray")
# }
# legend("bottomleft", legend = c("182 days", "273 days", "364 days", "448 days", "539 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
#
#
# ## Base STI higher-risk testing interval: 3003 - 10% annual screening
# ## Varying STI higher-risk testing interval: 3447:3451
# sims <- c(3447:3451, 3003)
# pal <- viridis::viridis(n = length(sims), option = "D")
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
#        mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
#        main = "HIV Incidence by Higher-Risk STI Screening Interval \n 50% HR / 10% Ann Cov",
#        xlab = "Week", ylab = "IR per 100 PYAR")
#   #abline(h = seq(0, 4, 0.2), lty = 2, col = "gray")
# }
#
# legend("bottomleft", legend = c("28 days", "91 days", "182 days", "273 days", "364 days", "Base - 10% Ann, 364 days"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()
#
# ## Histogram of partner numbers ------------------------------------------------
# tiff(filename = "analysis/Fig3b.tiff", height = 6, width = 11, units = "in", res = 250)
# #par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
#
# load("data/followup/sim.n3003.rda")
# zeropart <- mean(unname(colMeans(tail(sim$epi$zeropart, 26))))
# onepart <- mean(unname(colMeans(tail(sim$epi$onepart, 26))))
# twopart <- mean(unname(colMeans(tail(sim$epi$twopart, 26))))
# threepart <- mean(unname(colMeans(tail(sim$epi$threepart, 26))))
# fourpart <- mean(unname(colMeans(tail(sim$epi$fourpart, 26))))
# fivepart <- mean(unname(colMeans(tail(sim$epi$fivepart, 26))))
# sixpart <- mean(unname(colMeans(tail(sim$epi$sixpart, 26))))
# sevenpart <- mean(unname(colMeans(tail(sim$epi$sevenpart, 26))))
# eightpart <- mean(unname(colMeans(tail(sim$epi$eightpart, 26))))
# ninepart <- mean(unname(colMeans(tail(sim$epi$ninepart, 26))))
# tenpart <- mean(unname(colMeans(tail(sim$epi$tenpart, 26))))
# gttenpart <- mean(unname(colMeans(tail(sim$epi$gttenpart, 26))))
#
# Partners <- c(
#   rep(0, (100 * zeropart)),
#   rep(1, (100 * onepart)),
#   rep(2, (100 * twopart)),
#   rep(3, (100 * threepart)),
#   rep(4, (100 * fourpart)),
#   rep(5, (100 * fivepart)),
#   rep(6, (100 * sixpart)),
#   rep(7, (100 * sevenpart)),
#   rep(8, (100 * eightpart)),
#   rep(9, (100 * ninepart)),
#   rep(10, (100 * tenpart)),
#   rep(11, (100 * gttenpart))
#       )
#
# hist(Partners)
# dev.off()
#
