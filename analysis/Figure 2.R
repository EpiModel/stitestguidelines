## STI Testing Guidelines Figure 2 - Boxplots
##
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 10, type = "continuous")

## Coverage --------------------------------------------------------------
# 2 panels - HR ranging 0 to 100%, LR ranging 0 to 100%

# Baseline
load("data/followup/sim.n3000.rda")

sim.base <- sim
haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))

haz.sti <- as.numeric(colMeans(tail(sim.base$epi$ir100.sti, 52)))
ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
incid.base.sti <- unname(colSums(sim.base$epi$incid.sti[2:521, ]))

## Lower-risk coverage
# Sims 3001, 3003, 3005, 3007, 3009, 3011, 3013, 3015, 3017, 3019, 3021
# for lower-risk coverage - 0% higher-risk coverage
sims <- c(3001:3008)

df.hiv.pia <- data.frame(rep(NA, 256))
df.hivonly.nnt <- data.frame(rep(NA, 256))
df.hiv.nnt <- data.frame(rep(NA, 256))
df.gc.pia <- data.frame(rep(NA, 256))
df.gc.nnt <- data.frame(rep(NA, 256))
df.ct.pia <- data.frame(rep(NA, 256))
df.ct.nnt <- data.frame(rep(NA, 256))
df.syph.pia <- data.frame(rep(NA, 256))
df.syph.nnt <- data.frame(rep(NA, 256))
df.sti.pia <- data.frame(rep(NA, 256))
df.sti.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {

  load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))

  ir.comp <- unname(colMeans(sim$epi$ir100)) * 1000
  vec.hiv.nia <- round(ir.base - ir.comp, 1)
  df.hiv.pia[, i] <- vec.hiv.nia / ir.base

  ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  vec.gc.nia <- round(ir.base.gc - ir.comp.gc, 1)
  df.gc.pia[, i] <- vec.gc.nia / ir.base.gc

  ir.comp.ct <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  vec.ct.nia <- round(ir.base.ct - ir.comp.ct, 1)
  df.ct.pia[, i] <- vec.ct.nia / ir.base.ct

  ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
  vec.syph.nia <- round(ir.base.syph - ir.comp.syph, 1)
  df.syph.pia[, i] <- vec.syph.nia / ir.base.syph

  ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti)) * 1000
  vec.sti.nia <- round(ir.base.sti - ir.comp.sti, 1)
  df.sti.pia[, i] <- vec.sti.nia / ir.base.sti

  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (ir.base - ir.comp)
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (ir.base - ir.comp)
  df.gc.nnt[, i] <- (gc.asympt.tests) / (ir.base.gc - ir.comp.gc)
  df.ct.nnt[, i] <- (ct.asympt.tests) / (ir.base.ct - ir.comp.ct)
  df.syph.nnt[, i] <- (syph.asympt.tests) / (ir.base.syph - ir.comp.syph)
  df.sti.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (ir.base.sti - ir.comp.sti)

  cat("*")

}

names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  #c("20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")
  c("+5%", "+10%", "+15%", "+20%", "+25%", "+30%", "+35%", "+40%")

head(df.hiv.pia)
head(df.gc.pia)
head(df.ct.pia)
head(df.syph.pia)
head(df.sti.pia)
head(df.hiv.nnt)
head(df.gc.nnt)
head(df.ct.nnt)
head(df.syph.nnt)
head(df.sti.nnt)

## Higher-risk coverage
# Increases in higher-risk coverage at baseline coverage
# Baseline
load("data/followup/sim.n3000.rda")
sims <- c(3009, 3018, 3027, 3036, 3045, 3054, 3063,
          3072, 3081, 3090)

df2.hiv.pia <- data.frame(rep(NA, 256))
df2.hivonly.nnt <- data.frame(rep(NA, 256))
df2.hiv.nnt <- data.frame(rep(NA, 256))
df2.gc.pia <- data.frame(rep(NA, 256))
df2.gc.nnt <- data.frame(rep(NA, 256))
df2.ct.pia <- data.frame(rep(NA, 256))
df2.ct.nnt <- data.frame(rep(NA, 256))
df2.syph.pia <- data.frame(rep(NA, 256))
df2.syph.nnt <- data.frame(rep(NA, 256))
df2.sti.pia <- data.frame(rep(NA, 256))
df2.sti.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {

  load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))

  ir.comp <- unname(colMeans(sim$epi$ir100)) * 1000
  vec.hiv.nia <- round(ir.base - ir.comp, 1)
  df2.hiv.pia[, i] <- vec.hiv.nia / ir.base

  ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  vec.gc.nia <- round(ir.base.gc - ir.comp.gc, 1)
  df2.gc.pia[, i] <- vec.gc.nia / ir.base.gc

  ir.comp.ct <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  vec.ct.nia <- round(ir.base.ct - ir.comp.ct, 1)
  df2.ct.pia[, i] <- vec.ct.nia / ir.base.ct

  ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
  vec.syph.nia <- round(ir.base.syph - ir.comp.syph, 1)
  df2.syph.pia[, i] <- vec.syph.nia / ir.base.syph

  ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti)) * 1000
  vec.sti.nia <- round(ir.base.sti - ir.comp.sti, 1)
  df2.sti.pia[, i] <- vec.sti.nia / ir.base.sti

  hiv.tests <- unname(colSums(sim$epi$hivtests.nprep))
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests))
  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests))

  #HIV could be HIV tests or total STI tests
  df2.hivonly.nnt[, i] <- (hiv.tests) / (median(incid.base) - unname(colSums(sim$epi$incid)))
  df2.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (median(incid.base) - unname(colSums(sim$epi$incid)))
  df2.gc.nnt[, i] <- (gc.asympt.tests) / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
  df2.ct.nnt[, i] <- (ct.asympt.tests) / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))
  df2.syph.nnt[, i] <- (syph.asympt.tests) / (median(incid.base.syph) - unname(colSums(sim$epi$incid.syph)))
  df2.sti.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) /(median(incid.base.sti) - unname(colSums(sim$epi$incid.sti)))

  cat("*")

}

names(df2.hiv.pia) <- names(df2.gc.pia) <- names(df2.gc.nnt) <- names(df2.ct.pia) <- names(df2.sti.pia) <-
  names(df2.hivonly.nnt) <- names(df2.hiv.nnt) <- names(df2.ct.nnt) <- names(df2.syph.pia) <- names(df2.syph.nnt) <- names(df2.sti.nnt) <-
  c("10%", "20%", "30%", "40%", "50%", "60%", "70%",
    "80%", "90%", "100%")

head(df2.hiv.pia)
head(df2.gc.pia)
head(df2.ct.pia)
head(df2.syph.pia)
head(df2.sti.pia)
head(df2.hiv.nnt)
head(df2.gc.nnt)
head(df2.ct.nnt)
head(df2.syph.nnt)
head(df2.sti.nnt)


pal <- wes_palette("Zissou")[c(1, 5)]
tiff(filename = "analysis/Fig2.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))

# Combined STI
boxleft <- boxplot(df.sti.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(0, 1),
                   main = "PIA (STI) by coverage \n of sexually-active screening", las = 2,
                   xlab = "Coverage", ylab = "PIA",
                   cex.axis = 0.7)


boxright <- boxplot(df2.sti.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(0, 1),
                   main = "PIA (STI) by coverage \n of higher-risk screening", las = 2,
                   xlab = "Coverage", ylab = "PIA",
                   cex.axis = 0.7)

dev.off()

## Partner Cutoff --------------------------------------------------------------
#
# # Baseline
# load("data/followup/sim.n3003.rda")
# sim.base <- sim
# haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
# ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
# incid.base <- unname(colSums(sim.base$epi$incid))
#
# haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
# ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
# incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
#
# haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
# ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
# incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
#
# haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
# ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
# incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
#
# haz.sti <- as.numeric(colMeans(tail(sim.base$epi$ir100.sti, 52)))
# ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
# incid.base.sti <- unname(colSums(sim.base$epi$incid.sti[2:521, ]))
#
#
# # Sims 3494:3502 for Number of partners - 50% higher-risk coverage
# sims <- c(3494:3502)
# df.hiv.pia <- data.frame(rep(NA, 256))
# df.hivonly.nnt <- data.frame(rep(NA, 256))
# df.hiv.nnt <- data.frame(rep(NA, 256))
# df.gc.pia <- data.frame(rep(NA, 256))
# df.gc.nnt <- data.frame(rep(NA, 256))
# df.ct.pia <- data.frame(rep(NA, 256))
# df.ct.nnt <- data.frame(rep(NA, 256))
# df.syph.pia <- data.frame(rep(NA, 256))
# df.syph.nnt <- data.frame(rep(NA, 256))
# df.sti.pia <- data.frame(rep(NA, 256))
# df.sti.nnt <- data.frame(rep(NA, 256))
#
# for (i in seq_along(sims)) {
#
#   load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
#
#   ir.comp <- unname(colMeans(sim$epi$ir100)) * 1000
#   vec.hiv.nia <- round(ir.base - ir.comp, 1)
#   df.hiv.pia[, i] <- vec.hiv.nia / ir.base
#
#   ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
#   vec.gc.nia <- round(ir.base.gc - ir.comp.gc, 1)
#   df.gc.pia[, i] <- vec.gc.nia / ir.base.gc
#
#   ir.comp.ct <- unname(colMeans(sim$epi$ir100.gc)) * 1000
#   vec.ct.nia <- round(ir.base.ct - ir.comp.ct, 1)
#   df.ct.pia[, i] <- vec.ct.nia / ir.base.ct
#
#   ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
#   vec.syph.nia <- round(ir.base.syph - ir.comp.syph, 1)
#   df.syph.pia[, i] <- vec.syph.nia / ir.base.syph
#
#   ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti)) * 1000
#   vec.sti.nia <- round(ir.base.sti - ir.comp.sti, 1)
#   df.sti.pia[, i] <- vec.sti.nia / ir.base.sti
#
#   hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
#   gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
#   ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
#   syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))
#
#   #HIV could be HIV tests or total STI tests
#   df.hivonly.nnt[, i] <- (hiv.tests) / (ir.base - ir.comp)
#   df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (ir.base - ir.comp)
#   df.gc.nnt[, i] <- (gc.asympt.tests) / (ir.base.gc - ir.comp.gc)
#   df.ct.nnt[, i] <- (ct.asympt.tests) / (ir.base.ct - ir.comp.ct)
#   df.syph.nnt[, i] <- (syph.asympt.tests) / (ir.base.syph - ir.comp.syph)
#   df.sti.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (ir.base.sti - ir.comp.sti)
#
# }
# names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
#   names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
#   c("2", "3", "4", "5", "6", "7", "8", "9", "10")
#
# head(df.hiv.pia)
# head(df.gc.pia)
# head(df.ct.pia)
# head(df.syph.pia)
# head(df.sti.pia)
# head(df.hiv.nnt)
# head(df.gc.nnt)
# head(df.ct.nnt)
# head(df.syph.nnt)
# head(df.sti.nnt)
#
# pal <- wes_palette("Zissou")[c(1, 5)]
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 8, units = "in", res = 250)
# # HIV
# # tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# # par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# # boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
# #         main = "PIA (HIV) by Partner Cutoff", las = 2,
# #         xlab = "Partner Cutoff", ylab = "PIA",
# #         cex.axis = 0.7)
# #
# # boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (HIV) by Partner Cutoff", las = 2,
# #         xlab = "Partner Cutoff", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
# # dev.off()
#
# # GC
# # tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (NG) by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (NG) by Partner Cutoff", las = 2,
# #         xlab = "Partner Cutoff", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # CT
# # tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# # par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (CT) by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (CT) by Partner Cutoff",las = 2,
# #         xlab = "Partner Cutoff", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # Syph
# # tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# # par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (Syph) by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (Syph) by Partner Cutoff",las = 2,
# #         xlab = "Partner Cutoff", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # Combined STI
# boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (STI) by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "PIA",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()
#
#
# ## Screening Intervals - Lower Risk --------------------------------------------
# ## STI Testing Guidelines Figure 2 - Version 2 - Boxplots
#
# rm(list = ls())
# suppressMessages(library("EpiModelHIV"))
# library("EpiModelHPC")
# library("dplyr")
# source("analysis/fx.R")
#
# # Box Plots by Indications
# library('wesanderson')
# par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
# pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")
#
# # Baseline
# load("data/followup/sim.n3003.rda")
# sim.base <- sim
# incid.base <- unname(colSums(sim.base$epi$incid))
# incid.gc.base <- unname(colSums(sim.base$epi$incid.gc))
# incid.ct.base <- unname(colSums(sim.base$epi$incid.ct))
# incid.syph.base <- unname(colSums(sim.base$epi$incid.syph))
# #incid.sti.base <- unname(colSums(sim.base$epi$incid.sti))
# incid.sti.base <- incid.gc.base + incid.ct.base + incid.syph.base # temp workaround due to burnin having NA for time 5200
#
# sims <- c(3442:3446)
# df.hiv.pia <- data.frame(rep(NA, 256))
# df.hivonly.nnt <- data.frame(rep(NA, 256))
# df.hiv.nnt <- data.frame(rep(NA, 256))
# df.gc.pia <- data.frame(rep(NA, 256))
# df.gc.nnt <- data.frame(rep(NA, 256))
# df.ct.pia <- data.frame(rep(NA, 256))
# df.ct.nnt <- data.frame(rep(NA, 256))
# df.syph.pia <- data.frame(rep(NA, 256))
# df.syph.nnt <- data.frame(rep(NA, 256))
# df.sti.pia <- data.frame(rep(NA, 256))
# df.sti.nnt <- data.frame(rep(NA, 256))
#
# for (i in seq_along(sims)) {
#
#   load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
#
#   vec.hiv.nia <- round(incid.base - unname(colSums(sim$epi$incid)), 1)
#   df.hiv.pia[, i] <- vec.hiv.nia / incid.base
#
#   vec.gc.nia <- round(incid.gc.base - unname(colSums(sim$epi$incid.gc)), 1)
#   df.gc.pia[, i] <- vec.gc.nia / incid.gc.base
#
#   vec.ct.nia <- round(incid.ct.base - unname(colSums(sim$epi$incid.ct)), 1)
#   df.ct.pia[, i] <- vec.ct.nia / incid.ct.base
#
#   vec.syph.nia <- round(incid.syph.base - unname(colSums(sim$epi$incid.syph)), 1)
#   df.syph.pia[, i] <- vec.syph.nia / incid.syph.base
#
#   #vec.sti.nia <- round(incid.sti.base - unname(colSums(sim$epi$incid.sti)), 1)
#   vec.sti.nia <- round(incid.sti.base -
#                          (unname(colSums(sim$epi$incid.gc)) + unname(colSums(sim$epi$incid.ct)) +
#                             unname(colSums(sim$epi$incid.syph))), 1)
#   df.sti.pia[, i] <- vec.sti.nia / incid.sti.base
#
#
#   hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
#   gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
#   ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
#   syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))
#
#   #HIV could be HIV tests or total STI tests
#   df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
#   df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
#   df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
#   df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
#   df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
#   #df.sti.nnt[, i] <- syph.asympt.tests / (incid.sti.base - unname(colSums(sim$epi$incid.sti)))
#   df.sti.nnt[, i] <- syph.asympt.tests / (incid.sti.base -
#                                             (unname(colSums(sim$epi$incid.gc)) + unname(colSums(sim$epi$incid.ct)) +
#                                                unname(colSums(sim$epi$incid.syph))))
#
# }
# names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
#   names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
#   c("182 days", "273 days", "364 days", "448 days", "539 days")
#
# head(df.hiv.pia)
# head(df.gc.pia)
# head(df.ct.pia)
# head(df.syph.pia)
# head(df.sti.pia)
# head(df.hiv.nnt)
# head(df.gc.nnt)
# head(df.ct.nnt)
# head(df.syph.nnt)
# head(df.sti.nnt)
#
# pal <- wes_palette("Zissou")[c(1, 5)]
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
# # HIV
# # tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# # par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# # boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
# #         main = "PIA (HIV) by HR Interval", las = 2,
# #         xlab = "Interval", ylab = "PIA",
# #         cex.axis = 0.7)
# #
# # boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (HIV) by HR Interval", las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% Ann Cov (6 months), 0% HR Cov", outer = TRUE)
# # dev.off()
#
# # GC
# # tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (NG) by HR Interval", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (NG) by HR Interval", las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% Ann Cov (6 months), 0% HR Cov", outer = TRUE)
# # dev.off()
#
# # CT
# # tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# # par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (CT) by HR Interval", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (CT) by HR Interval",las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% Ann Cov (6 months), 0% HR Cov", outer = TRUE)
# # dev.off()
#
# # Syph
# # tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# # par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (Syph) by HR Interval", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (Syph) by HR Interval",las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% Ann Cov (6 months), 0% HR Cov", outer = TRUE)
# # dev.off()
#
# # Combined STI
# boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (STI) by Partner Cutoff", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
# #         main = "NNT (STI) by Partner Cutoff", las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# title("50% Ann Cov (6 months), 0% HR Cov", outer = TRUE)
# dev.off()
#
#
#
# ## Screening Intervals - Higher-Risk -------------------------------------------
# ## STI Testing Guidelines Figure 2 - Version 2 - Boxplots
#
# rm(list = ls())
# suppressMessages(library("EpiModelHIV"))
# library("EpiModelHPC")
# library("dplyr")
# source("analysis/fx.R")
#
# # Box Plots by Indications
# library('wesanderson')
# par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
# pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")
#
# # Baseline
# load("data/followup/sim.n3003.rda")
# sim.base <- sim
# incid.base <- unname(colSums(sim.base$epi$incid))
# incid.gc.base <- unname(colSums(sim.base$epi$incid.gc))
# incid.ct.base <- unname(colSums(sim.base$epi$incid.ct))
# incid.syph.base <- unname(colSums(sim.base$epi$incid.syph))
# #incid.sti.base <- unname(colSums(sim.base$epi$incid.sti))
# incid.sti.base <- incid.gc.base + incid.ct.base + incid.syph.base # temp workaround due to burnin having NA for time 5200
#
#
# sims <- c(3447:3451)
# df.hiv.pia <- data.frame(rep(NA, 256))
# df.hivonly.nnt <- data.frame(rep(NA, 256))
# df.hiv.nnt <- data.frame(rep(NA, 256))
# df.gc.pia <- data.frame(rep(NA, 256))
# df.gc.nnt <- data.frame(rep(NA, 256))
# df.ct.pia <- data.frame(rep(NA, 256))
# df.ct.nnt <- data.frame(rep(NA, 256))
# df.syph.pia <- data.frame(rep(NA, 256))
# df.syph.nnt <- data.frame(rep(NA, 256))
# df.sti.pia <- data.frame(rep(NA, 256))
# df.sti.nnt <- data.frame(rep(NA, 256))
#
# for (i in seq_along(sims)) {
#
#   sim <- truncate_sim(sim.base, at = 5201)
#
#   load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
#
#   vec.hiv.nia <- round(incid.base - unname(colSums(sim$epi$incid)), 1)
#   df.hiv.pia[, i] <- vec.hiv.nia / incid.base
#
#   vec.gc.nia <- round(incid.gc.base - unname(colSums(sim$epi$incid.gc)), 1)
#   df.gc.pia[, i] <- vec.gc.nia / incid.gc.base
#
#   vec.ct.nia <- round(incid.ct.base - unname(colSums(sim$epi$incid.ct)), 1)
#   df.ct.pia[, i] <- vec.ct.nia / incid.ct.base
#
#   vec.syph.nia <- round(incid.syph.base - unname(colSums(sim$epi$incid.syph)), 1)
#   df.syph.pia[, i] <- vec.syph.nia / incid.syph.base
#
#   #vec.sti.nia <- round(incid.sti.base - unname(colSums(sim$epi$incid.sti)), 1)
#   vec.sti.nia <- round(incid.sti.base -
#                          (unname(colSums(sim$epi$incid.gc)) + unname(colSums(sim$epi$incid.ct)) +
#                             unname(colSums(sim$epi$incid.syph))), 1)
#   df.sti.pia[, i] <- vec.sti.nia / incid.sti.base
#
#   hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
#   gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
#   ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
#   syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))
#
#   #HIV could be HIV tests or total STI tests
#   df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
#   df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
#   df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
#   df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
#   df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
#   #df.sti.nnt[, i] <- syph.asympt.tests / (incid.sti.base - unname(colSums(sim$epi$incid.sti)))
#   df.sti.nnt[, i] <- syph.asympt.tests / (incid.sti.base -
#                                             (unname(colSums(sim$epi$incid.gc)) + unname(colSums(sim$epi$incid.ct)) +
#                                                unname(colSums(sim$epi$incid.syph))))
#
# }
# names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
#   names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
#   c("28 days", "91 days", "182 days", "237 days", "364 days")
#
# head(df.hiv.pia)
# head(df.gc.pia)
# head(df.ct.pia)
# head(df.syph.pia)
# head(df.sti.pia)
# head(df.hiv.nnt)
# head(df.gc.nnt)
# head(df.ct.nnt)
# head(df.syph.nnt)
# head(df.sti.nnt)
#
# pal <- wes_palette("Zissou")[c(1, 5)]
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# # HIV
# # tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# # par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# # boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
# #         main = "PIA (HIV) by HR Interval", las = 2,
# #         xlab = "Interval", ylab = "PIA",
# #         cex.axis = 0.7)
# #
# # boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (HIV) by HR Interval", las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # GC
# # tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (NG) by HR Interval", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (NG) by HR Interval", las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # CT
# # tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# # par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (CT) by HR Interval", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (CT) by HR Interval",las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # Syph
# # tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# # par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (Syph) by HR Interval", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)),
# #         main = "NNT (Syph) by HR Interval",las = 2,
# #         xlab = "Interval", ylab = "NNT"
# #         cex.axis = 0.7)
# # title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# # dev.off()
#
# # Combined STI
# boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
#         main = "PIA (STI) by Partner Cutoff", las = 2,
#         xlab = "Interval", ylab = "PIA",
#         cex.axis = 0.7)
# # boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
# #         col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
# #         main = "NNT (STI) by Partner Cutoff", las = 2,
# #         xlab = "Interval", ylab = "NNT",
# #         cex.axis = 0.7)
# #
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()
