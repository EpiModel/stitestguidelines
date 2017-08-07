#
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
