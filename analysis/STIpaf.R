## STI PAF Figures

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Base STI relative risks: 3000
## Base no relative risks: 2999

#tiff(filename = "analysis/Fig3b.tiff", height = 6, width = 11, units = "in", res = 250)
tiff(filename = "analysis/Fig4.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(2999, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100", add = i > 1,
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3,
         main = "HIV incidence by STI Relative Risk",
         xlab = "Week", ylab = "IR per 100 PYAR")
}
legend("bottomleft", legend = c("No increased risk of acquisition", "Current Relative Risk"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()

par(mfrow = c(1,1), oma = c(0,0,2,0))
plot(sim, y = "stiactiveind", mean.col = "purple", ylim = c(0, 1))#, xlim = c(0, sim$control$nsteps - sim$param$stitest.start), ylim = c(0, 1))
plot(sim, y = "newpartner", mean.col = "orange", add = TRUE)
plot(sim, y = "recentpartners", mean.col = "green", add = TRUE)
plot(sim, y = "concurrpart", mean.col = "blue", add = TRUE)
plot(sim, y = "partnersti", mean.col = "brown", add = TRUE)
plot(sim, y = "uai.nmain", mean.col = "black", add = TRUE)
plot(sim, y = "uai.any", mean.col = "gray", add = TRUE)
plot(sim, y = "recentSTI", mean.col = "red", add = TRUE)
abline(h = c(seq(0.1, 0.9, 0.1)), lty = 2, col = "gray")
legend(400, 0.9, lty = c(rep(1, 8)), 
       col = c("purple", "orange", "green", "blue", "brown", "black", "gray", "red"),
       c("Sexually Active", "New Partner", ">1 Recent Partners", "Partner is Concurrent", 
         "Partner had STI", "CAI in Non-main", "Any CAI", "Recent STI"))
title("STI Testing Indications")

par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
plot(sim, y = "hiv_sum", mean.col = "blue", qnts.col = "blue", qnts.alpha = 0.2)
plot(sim, y = "sti_hiv_sum", mean.col = "green", qnts.col = "green", qnts.alpha = 0.2, add = TRUE)
plot(sim, y = "sti_u_hiv_sum", mean.col = "red", qnts.col = "red", qnts.alpha = 0.2,add = TRUE)
plot(sim, y = "sti_r_hiv_sum", mean.col = "orange", qnts.col = "orange", qnts.alpha = 0.2,add = TRUE)
plot(sim, y = "sti_syph_hiv_sum", mean.col = "purple", qnts.col = "purple", qnts.alpha = 0.2,add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("blue", "green", "red", "orange", "purple"), 
       c("HIV Sum", "STI & HIV sum", "USTI & HIV sum", "RSTI & HIV sum", "Syph & HIV sum"))
title("HIV Infections per time step")

par(mfrow = c(1, 1))
plot(sim, y = "sti_paf", mean.col = "blue", qnts.col = "blue", qnts = c(0.5), qnts.alpha = 0.05)
abline(h = c(seq(0.1, 0.9, 0.1)), lty = 2, col = "gray")
plot(sim, y = "sti_u_paf", mean.col = "green", qnts.col = "green", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sti_r_paf", mean.col = "red", qnts.col = "red", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sti_syph_paf", mean.col = "orange", qnts.col = "orange", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1), col = c("blue", "green", "red", "orange"), 
       c("Total STI PAF", "Urethral NG/CT PAF", "Rectal NG/CT PAF", "Syph PAF"))
title("STI Population Attributable Fraction (PAF) and IQR for HIV Infection")

plot(sim, y = "sti_u_paf", mean.col = "blue")
plot(sim, y = "sti_u_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_u_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("U STI PAF", "U Sympt STI PAF", "U Asympt STI PAF"))

plot(sim, y = "sti_r_paf", mean.col = "blue")
plot(sim, y = "sti_r_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_r_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("R STI PAF", "R Sympt STI PAF", "R Asympt STI PAF"))

plot(sim, y = "sti_syph_paf", mean.col = "blue")
plot(sim, y = "sti_syph_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_syph_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("Syph PAF", "Syph Sympt PAF", "Syph Asympt PAF"))