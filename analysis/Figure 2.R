## STI Testing Guidelines Figure 2 - Boxplots
## Partner Cutoff --------------------------------------------------------------

rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")

# Baseline
load("data/followup/sim.n3003.rda")
sim.base <- sim
incid.base <- unname(colSums(sim.base$epi$incid))
incid.gc.base <- unname(colSums(sim.base$epi$incid.gc))
incid.ct.base <- unname(colSums(sim.base$epi$incid.ct))
incid.syph.base <- unname(colSums(sim.base$epi$incid.syph))
incid.sti.base <- unname(colSums(sim.base$epi$incid.sti))

# Sims 3494:3502 for Number of partners - 50% higher-risk coverage
sims <- c(3003, 3494:3502)
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

  vec.hiv.nia <- round(incid.base - unname(colSums(sim$epi$incid)), 1)
  df.hiv.pia[, i] <- vec.hiv.nia / incid.base

  vec.gc.nia <- round(incid.gc.base - unname(colSums(sim$epi$incid.gc)), 1)
  df.gc.pia[, i] <- vec.gc.nia / incid.gc.base

  vec.ct.nia <- round(incid.ct.base - unname(colSums(sim$epi$incid.ct)), 1)
  df.ct.pia[, i] <- vec.ct.nia / incid.ct.base

  vec.syph.nia <- round(incid.syph.base - unname(colSums(sim$epi$incid.syph)), 1)
  df.syph.pia[, i] <- vec.syph.nia / incid.syph.base

  vec.sti.nia <- round(incid.sti.base - unname(colSums(sim$epi$incid.sti)), 1)
  df.sti.pia[, i] <- vec.sti.nia / incid.sti.base

  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))

}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

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

pal <- wes_palette("Zissou")[c(1, 5)]
# HIV
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
#         main = "PIA by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "Proportion HIV Infections Averted",
#         cex.axis = 0.7)
#
# boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "HIV Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
# dev.off()

# GC
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion NG Infections Averted",
        cex.axis = 0.7)
# boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "NG Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# CT
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion CT Infections Averted",
        cex.axis = 0.7)
# boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "CT Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Syph
# tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion Syph Infections Averted",
        cex.axis = 0.7)
# boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "Syph Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Combined STI
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
main = "PIA by Partner Cutoff", las = 2,
xlab = "Partner Cutoff", ylab = "Proportion STI Infections Averted",
cex.axis = 0.7)
title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
dev.off()


## Screening Intervals - Lower Risk --------------------------------------------
## STI Testing Guidelines Figure 2 - Version 2 - Boxplots

rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")

# Baseline
load("data/followup/sim.n3003.rda")
sim.base <- sim
incid.base <- unname(colSums(sim.base$epi$incid))
incid.gc.base <- unname(colSums(sim.base$epi$incid.gc))
incid.ct.base <- unname(colSums(sim.base$epi$incid.ct))
incid.syph.base <- unname(colSums(sim.base$epi$incid.syph))
incid.sti.base <- unname(colSums(sim.base$epi$incid.sti))

sims <- c(3442:3446)
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

  vec.hiv.nia <- round(incid.base - unname(colSums(sim$epi$incid)), 1)
  df.hiv.pia[, i] <- vec.hiv.nia / incid.base

  vec.gc.nia <- round(incid.gc.base - unname(colSums(sim$epi$incid.gc)), 1)
  df.gc.pia[, i] <- vec.gc.nia / incid.gc.base

  vec.ct.nia <- round(incid.ct.base - unname(colSums(sim$epi$incid.ct)), 1)
  df.ct.pia[, i] <- vec.ct.nia / incid.ct.base

  vec.syph.nia <- round(incid.syph.base - unname(colSums(sim$epi$incid.syph)), 1)
  df.syph.pia[, i] <- vec.syph.nia / incid.syph.base

  vec.sti.nia <- round(incid.sti.base - unname(colSums(sim$epi$incid.sti)), 1)
  df.sti.pia[, i] <- vec.sti.nia / incid.sti.base


  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))

}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

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

pal <- wes_palette("Zissou")[c(1, 5)]
# HIV
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
#         main = "PIA by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "Proportion HIV Infections Averted",
#         cex.axis = 0.7)
#
# boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "HIV Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# GC
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion NG Infections Averted",
        cex.axis = 0.7)
# boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "NG Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# CT
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion CT Infections Averted",
        cex.axis = 0.7)
# boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "CT Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Syph
# tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion Syph Infections Averted",
        cex.axis = 0.7)
# boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "Syph Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Combined STI
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion STI Infections Averted",
        cex.axis = 0.7)
title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
dev.off()



## Screening Intervals - Higher-Risk -------------------------------------------
## STI Testing Guidelines Figure 2 - Version 2 - Boxplots

rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")

# Baseline
load("data/followup/sim.n3003.rda")
sim.base <- sim
incid.base <- unname(colSums(sim.base$epi$incid))
incid.gc.base <- unname(colSums(sim.base$epi$incid.gc))
incid.ct.base <- unname(colSums(sim.base$epi$incid.ct))
incid.syph.base <- unname(colSums(sim.base$epi$incid.syph))
incid.sti.base <- unname(colSums(sim.base$epi$incid.sti))

sims <- c(3447:3451)
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

  vec.hiv.nia <- round(incid.base - unname(colSums(sim$epi$incid)), 1)
  df.hiv.pia[, i] <- vec.hiv.nia / incid.base

  vec.gc.nia <- round(incid.gc.base - unname(colSums(sim$epi$incid.gc)), 1)
  df.gc.pia[, i] <- vec.gc.nia / incid.gc.base

  vec.ct.nia <- round(incid.ct.base - unname(colSums(sim$epi$incid.ct)), 1)
  df.ct.pia[, i] <- vec.ct.nia / incid.ct.base

  vec.syph.nia <- round(incid.syph.base - unname(colSums(sim$epi$incid.syph)), 1)
  df.syph.pia[, i] <- vec.syph.nia / incid.syph.base

  vec.sti.nia <- round(incid.sti.base - unname(colSums(sim$epi$incid.sti)), 1)
  df.sti.pia[, i] <- vec.sti.nia / incid.sti.base

  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))

}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

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

pal <- wes_palette("Zissou")[c(1, 5)]
# HIV
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
#         main = "PIA by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "Proportion HIV Infections Averted",
#         cex.axis = 0.7)
#
# boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "HIV Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
# dev.off()

# GC
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion NG Infections Averted",
        cex.axis = 0.7)
# boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "NG Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# CT
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion CT Infections Averted",
        cex.axis = 0.7)
# boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "CT Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Syph
# tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion Syph Infections Averted",
        cex.axis = 0.7)
# boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "Syph Number Needed to Treat",
#         cex.axis = 0.7)
# title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Combined STI
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion STI Infections Averted",
        cex.axis = 0.7)
title("50% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
dev.off()
