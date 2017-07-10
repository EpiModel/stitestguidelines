## STI PAF Plot

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
library("wesanderson")
library("viridis")
source("analysis/fx.R")

tiff(filename = "analysis/HIVTrans1.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1, 1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3000, 4001:4008,4010)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    par(mfrow = c(1, 1), oma = c(3, 0, 2, 0))
    plot(sim, y = "ir100", add = i > 1, ylim = c(0, 8),
         mean.col = pal[i], qnts = 0.5, qnts.col = pal[i], qnts.alpha = 0.1,
         main = "HIV incidence and IQR by STI Relative Transmission Risk",
         xlab = "Week", ylab = "IR per 100 PYAR")
}
legend(title = "Relative Risks", "bottomleft",
       legend = c("NG/CT = 1.0, Syph = 1.0", "NG/CT = 1.1, Syph = 1.0", "NG/CT = 1.2, Syph = 1.0",
                  "NG/CT = 1.3, Syph = 1.0","NG/CT = 1.4, Syph = 1.0", "NG/CT = 1.5, Syph = 1.0",
                  "NG/CT = 1.6, Syph = 1.0", "NG/CT = 1.7, Syph = 1.0","NG/CT = 1.8, Syph = 1.0",
                  "NG/CT = 2.0, Syph = 1.0"),
       ncol = 2, col = pal, lwd = 3, cex = 0.85, bty = "n")
mtext(side = 1, text = "Relative Risks for acquisition: Syphilis (Receptive) = 1.8029, Syphilis (Insertive) = 1.1989,
      Rectal Gonorrhea = 1.8029, Urethral Gonorrhea = 1.1989, Rectal Chlamydia = 1.8029, Urethral Chlamydia = 1.1989",
      at = 0.5, padj = 1, outer = TRUE)

dev.off()

# tiff(filename = "analysis/FigPAF2.tiff", height = 6, width = 11, units = "in", res = 250)
# par(mfrow = c(1, 1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
# sims <- c(3000, 4001:4008, 4010)
# pal <- viridis::viridis(n = length(sims), option = "D")
#
# for (i in seq_along(sims)) {
#   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
#   load(fn)
#   plot(sim, y = "ir100.sti", add = i > 1, mean.smooth = FALSE,
#        mean.col = pal[i], qnts = 0.5, qnts.col = pal[i], qnts.alpha = 0.1,
#        main = "Combined STI incidence by STI Relative Transmission Risk",
#        xlab = "Week", ylab = "IR per 100 PYAR")
# }
# legend("bottomleft", legend = c("NG/CT = 1", "NG/CT = 1.1", "NG/CT = 1.2",
#                                 "NG/CT = 1.3","NG/CT = 1.4", "NG/CT = 1.5",
#                                 "NG/CT = 1.6", "NG/CT = 1.7","NG/CT = 1.8",
#                                 "NG/CT = 2.0"),
#        col = pal, lwd = 3, cex = 0.85, bty = "n")
# dev.off()

tiff(filename = "analysis/HIVTrans2.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1, 1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3000, 4021, 4032, 4043, 4054, 4065, 4076, 4087, 4098, 4109, 4120)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  par(mfrow = c(1, 1), oma = c(3, 0, 2, 0))
  plot(sim, y = "ir100", add = i > 1, ylim = c(0, 8),
       mean.col = pal[i], qnts = 0.5, qnts.col = pal[i], qnts.alpha = 0.1,
       main = "HIV incidence and IQR by STI Relative Transmission Risk",
       xlab = "Week", ylab = "IR per 100 PYAR")
}
legend(title = "Relative Risks", "bottomleft",
       legend = c("Syph/NG/CT = 1.0", "Syph/NG/CT = 1.1", "Syph/NG/CT = 1.2",
                  "Syph/NG/CT = 1.3", "Syph/NG/CT = 1.4", "Syph/NG/CT = 1.5",
                  "Syph/NG/CT = 1.6", "Syph/NG/CT = 1.7", "Syph/NG/CT = 1.8",
                  "Syph/NG/CT = 1.9", "Syph/NG/CT = 2.0"),
       ncol = 2, col = pal, lwd = 3, cex = 0.85, bty = "n")
mtext(side = 1, text = "Relative Risks for acquisition: Syphilis (Receptive) = 1.8029, Syphilis (Insertive) = 1.1989,
      Rectal Gonorrhea = 1.8029, Urethral Gonorrhea = 1.1989, Rectal Chlamydia = 1.8029, Urethral Chlamydia = 1.1989",
      at = 0.5, padj = 1, outer = TRUE)

dev.off()

tiff(filename = "analysis/HIVTransprev.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1, 1), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3000, 4021, 4032, 4043, 4054, 4065, 4076, 4087, 4098, 4109, 4120)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  par(mfrow = c(1, 1), oma = c(3, 0, 2, 0))
  plot(sim, y = "i.prev", add = i > 1, ylim = c(0, 0.3),
       mean.col = pal[i], qnts = 0.5, qnts.col = pal[i], qnts.alpha = 0.1,
       main = "HIV Prevalence and IQR by STI Relative Transmission Risk",
       xlab = "Week", ylab = "Prevalence")
}
legend(title = "Relative Risks", "bottomleft",
       legend = c("Syph/NG/CT = 1.0", "Syph/NG/CT = 1.1", "Syph/NG/CT = 1.2",
                  "Syph/NG/CT = 1.3", "Syph/NG/CT = 1.4", "Syph/NG/CT = 1.5",
                  "Syph/NG/CT = 1.6", "Syph/NG/CT = 1.7", "Syph/NG/CT = 1.8",
                  "Syph/NG/CT = 1.9", "Syph/NG/CT = 2.0"),
       ncol = 2, col = pal, lwd = 3, cex = 0.85, bty = "n")
mtext(side = 1, text = "Relative Risks for acquisition: Syphilis (Receptive) = 1.8029, Syphilis (Insertive) = 1.1989,
      Rectal Gonorrhea = 1.8029, Urethral Gonorrhea = 1.1989, Rectal Chlamydia = 1.8029, Urethral Chlamydia = 1.1989",
      at = 0.5, padj = 1, outer = TRUE)

dev.off()


######### Prevalence of STIs by serostatus
# GC
tiff(filename = "analysis/GCbyHIV.tiff", height = 6, width = 11, units = "in", res = 250)
rm(list = ls())
load("data/followup/sim.n3000.rda")
par(mfrow = c(2, 2), oma = c(0,0,0,0))
plot(sim, y = "prev.gc.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.gc.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("NG and HIV-positive", "NG and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title("NG Prevalence by HIV status")

plot(sim, y = "prev.rgc.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.rgc.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", c(legend = "Rectal NG and HIV-positive", "Rectal NG and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title("")
title("Rectal NG Prevalence by HIV status")

plot(sim, y = "prev.ugc.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.ugc.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("Urethral NG and HIV-positive", "Urethral NG and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title("")
title("Urethral NG Prevalence by HIV status")
dev.off()

# CT
tiff(filename = "analysis/CTbyHIV.tiff", height = 6, width = 11, units = "in", res = 250)
rm(list = ls())
load("data/followup/sim.n3000.rda")
par(mfrow = c(2, 2), oma = c(0,0,2,0))
plot(sim, y = "prev.ct.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.ct.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("CT and HIV-positive", "CT and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title()
title("CT Prevalence by HIV status")

plot(sim, y = "prev.rct.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.rct.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("Rectal CT and HIV-positive", "Rectal CT and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title()
title("Rectal CT Prevalence by HIV status")

plot(sim, y = "prev.uct.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.uct.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("Urethral CT and HIV-positive", "Urethral CT and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title()
title("Urethral CT Prevalence by HIV status")
dev.off()

# Syphilis
tiff(filename = "analysis/SyphbyHIV.tiff", height = 6, width = 11, units = "in", res = 250)
rm(list = ls())
load("data/followup/sim.n3000.rda")
par(mfrow = c(1, 2), oma = c(0,0,2,0))
plot(sim, y = "prev.syph.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.syph.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("Syph and HIV-positive", "Syph and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title()
title("Syphilis Prevalence by HIV status")

plot(sim, y = "prev.primsecosyph.hivpos", ylab = "Prevalence", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.primsecosyph.hivneg", mean.col = "red", qnts.col = "red", add = TRUE)
legend("topleft", legend = c("P&S Syph and HIV-positive", "P&S Syph and HIV-negative"), lty = c(1, 1), col = c("blue", "red"))
title()
title("Primary and Secondary Syphilis Prevalence by HIV status")
dev.off()

####### PAF graph
tiff(filename = "analysis/PAF.tiff", height = 6, width = 11, units = "in", res = 250)
rm(list = ls())
load("data/followup/sim.n3000.rda")
par(mfrow = c(1, 1), oma = c(3, 0, 2, 0))
plot(sim, y = "sum_GC", mean.col = "blue", qnts.col = "blue", qnts = c(0.5), qnts.alpha = 0.05)
abline(h = c(seq(0.1, 0.9, 0.1)), lty = 2, col = "gray")
plot(sim, y = "sum_CT", mean.col = "green", qnts.col = "green", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sum_syph", mean.col = "red", qnts.col = "red", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sum_rectal", mean.col = "orange", qnts.col = "orange", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sum_urethral", mean.col = "purple", qnts.col = "purple", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sti_syph_paf", mean.col = "orange", qnts.col = "orange", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1), col = c("blue", "green", "red", "orange", "purple"),
       legend = c("NG", "CT", "Syphilis", "Urethral NG/CT", "Rectal NG/CT"))
title()
title("Prevalence of Site-Specific STI among Incident HIV Infections")
mtext(side = 1, text = "Relative Risks for acquisition: Syphilis (Receptive) = 1.8029, Syphilis (Insertive) = 1.1989,
      Rectal Gonorrhea = 1.8029, Urethral Gonorrhea = 1.1989, Rectal Chlamydia = 1.8029, Urethral Chlamydia = 1.1989",
      at = 0.5, padj = 1, outer = TRUE)
dev.off()

# plot(sim, y = "sti_u_paf", mean.col = "blue")
# plot(sim, y = "sti_u_sympt_paf", mean.col = "green", add = TRUE)
# plot(sim, y = "sti_u_asympt_paf", mean.col = "red", add = TRUE)
# legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"),
#        c("U STI PAF", "U Sympt STI PAF", "U Asympt STI PAF"))
#
# plot(sim, y = "sti_r_paf", mean.col = "blue")
# plot(sim, y = "sti_r_sympt_paf", mean.col = "green", add = TRUE)
# plot(sim, y = "sti_r_asympt_paf", mean.col = "red", add = TRUE)
# legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"),
#        c("R STI PAF", "R Sympt STI PAF", "R Asympt STI PAF"))
#
# plot(sim, y = "sti_syph_paf", mean.col = "blue")
# plot(sim, y = "sti_syph_sympt_paf", mean.col = "green", add = TRUE)
# plot(sim, y = "sti_syph_asympt_paf", mean.col = "red", add = TRUE)
# legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"),
#        c("Syph PAF", "Syph Sympt PAF", "Syph Asympt PAF"))







