## STI Testing Guidelines Figure 1 - Boxplots

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
load("data/sim.n3000.rda")
#load("data/followup/sim.n3000.rda")

sim.base <- sim
haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))
tests.base <- unname(colSums(sim.base$epi$hivtests.nprep))

haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests))

haz.sti <- as.numeric(colMeans(tail(sim.base$epi$ir100.sti, 52)))
ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
incid.base.sti <- unname(colSums(sim.base$epi$incid.sti[2:521, ]))
tests.sti.base <- unname(colSums(sim.base$epi$stiasympttests))

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

  #load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
  load(list.files("data/", pattern = as.character(sims[i]), full.names = TRUE))

  incid <- unname(colSums(sim$epi$incid))
  vec.hiv.nia <- round(incid.base - incid, 1)
  df.hiv.pia[, i] <- vec.hiv.nia / ir.base

  incid.gc <- unname(colSums(sim$epi$incid.gc))
  vec.gc.nia <- incid.base.gc - incid.gc
  df.gc.pia[, i] <- vec.gc.nia / incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.ct))
  vec.ct.nia <- incid.base.ct - incid.ct
  df.ct.pia[, i] <- vec.ct.nia / incid.base.ct

  incid.syph <- unname(colSums(sim$epi$incid.syph))
  vec.syph.nia <- incid.base.syph - incid.syph
  df.syph.pia[, i] <- vec.syph.nia / incid.base.syph

  incid.sti <- unname(colSums(sim$epi$incid.sti))
  vec.sti.nia <- incid.base.sti - incid.sti
  df.sti.pia[, i] <- vec.sti.nia / incid.base.sti

  # Tests
  hiv.tests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
  sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))

  #HIV could be HIV tests or total STI tests
  df.hiv.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- (sti.asympt.tests  - tests.sti.base) / (incid.base.sti - unname(colSums(sim$epi$incid.sti)))

  cat("*")

}

names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
 names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
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
sims <- c(3018, 3036, 3054, 3072, 3090, 3108, 3126, 3144, 3162, 3180)

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


  incid <- unname(colSums(sim$epi$incid))
  vec.hiv.nia <- round(incid.base - incid, 1)
  df2.hiv.pia[, i] <- vec.hiv.nia / incid.base

  incid.gc <- unname(colSums(sim$epi$incid.gc))
  vec.gc.nia <- incid.base.gc - incid.gc
  df2.gc.pia[, i] <- vec.gc.nia / incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.ct))
  vec.ct.nia <- incid.base.ct - incid.ct
  df2.ct.pia[, i] <- vec.ct.nia / incid.base.ct

  incid.syph <- unname(colSums(sim$epi$incid.syph))
  vec.syph.nia <- incid.base.syph - incid.syph
  df2.syph.pia[, i] <- vec.syph.nia / incid.base.syph

  incid.sti <- unname(colSums(sim$epi$incid.sti))
  vec.sti.nia <- incid.base.sti - incid.sti
  df2.sti.pia[, i] <- vec.sti.nia / incid.base.sti

  # Tests
  hiv.tests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
  sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))

  #HIV could be HIV tests or total STI tests
  df2.hiv.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df2.gc.nnt[, i] <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - unname(colSums(sim$epi$incid.gc)))
  df2.ct.nnt[, i] <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - unname(colSums(sim$epi$incid.ct)))
  df2.syph.nnt[, i] <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - unname(colSums(sim$epi$incid.syph)))
  df2.sti.nnt[, i] <- (sti.asympt.tests  - tests.sti.base) / (incid.base.sti - unname(colSums(sim$epi$incid.sti)))

  cat("*")

  cat("*")

}

names(df2.hiv.pia) <- names(df2.gc.pia) <- names(df2.gc.nnt) <- names(df2.ct.pia) <- names(df2.sti.pia) <-
  names(df2.hiv.nnt) <- names(df2.ct.nnt) <- names(df2.syph.pia) <- names(df2.syph.nnt) <- names(df2.sti.nnt) <-
  c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%")

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

pal <- wes_palette("GrandBudapest2")[2]

# Combined STI -----------------------------------------------------------------
tiff(filename = "analysis/Fig1PIA.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))

# Combined STI
boxleft <- boxplot(df.sti.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                   main = "Percent of total STI infections averted (PIA) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Percent of Infections Averted (PIA)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.sti.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                   main = "Percent of total STI infections averted (PIA) \n by coverage of higher-risk screening", las = 2,
                   xlab = "Coverage of higher-risk screening", ylab = "Percent of Infections Averted (PIA)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()

tiff(filename = "analysis/Fig1NNT.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))

boxleft <- boxplot(df.sti.nnt, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)),
                   main = "Number Needed to Treat (NNT) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Number Needed to Treat (NNT)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.sti.nnt, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)),
                    main = "Number Needed to Treat (NNT) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Number Needed to Treat (NNT)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()


# NG  -----------------------------------------------------------------
tiff(filename = "analysis/Fig1PIAGC.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))


boxleft <- boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                   main = "Percent of total NG infections averted (PIA) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Percent of Infections Averted (PIA)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.gc.pia, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                    main = "Percent of total NG infections averted (PIA) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Percent of Infections Averted (PIA)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()

tiff(filename = "analysis/Fig1NNTGC.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))
boxleft <- boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)),
                   main = "NG Number Needed to Treat (NNT) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Number Needed to Treat (NNT)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.gc.nnt, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)),
                    main = "NG Number Needed to Treat (NNT) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Number Needed to Treat (NNT)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()

# CT  -----------------------------------------------------------------
tiff(filename = "analysis/Fig1PIACT.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))


boxleft <- boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                   main = "Percent of total CT infections averted (PIA) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Percent of Infections Averted (PIA)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.ct.pia, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                    main = "Percent of total CT infections averted (PIA) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Percent of Infections Averted (PIA)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()

tiff(filename = "analysis/Fig1NNTCT.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))
boxleft <- boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)),
                   main = "CT Number Needed to Treat (NNT) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Number Needed to Treat (NNT)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.ct.nnt, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)),
                    main = "CT Number Needed to Treat (NNT) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Number Needed to Treat (NNT)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()

# Syph  -----------------------------------------------------------------
tiff(filename = "analysis/Fig1PIASyph.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))


boxleft <- boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                   main = "Percent of total Syph infections averted (PIA) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Percent of Infections Averted (PIA)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.syph.pia, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)), ylim = c(-1, 1),
                    main = "Percent of total Syph infections averted (PIA) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Percent of Infections Averted (PIA)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()

tiff(filename = "analysis/Fig1NNTSyph.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(1, 1, 2, 1), mgp = c(3, 0.75, 0))

boxleft <- boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
                   col = c(rep(pal[1], 10)),
                   main = "Syph Number Needed to Treat (NNT) \n by coverage of sexually-active screening", las = 2,
                   xlab = "Coverage of sexually active screening", ylab = "Number Needed to Treat (NNT)",
                   cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

boxright <- boxplot(df2.syph.nnt, outline = FALSE, medlwd = 1.1,
                    col = c(rep(pal[1], 10)),
                    main = "Syph Number Needed to Treat (NNT) \n by coverage of higher-risk screening", las = 2,
                    xlab = "Coverage of higher-risk screening", ylab = "Number Needed to Treat (NNT)",
                    cex.axis = 0.7, cex.main = 0.8, cex.lab = 0.8)

dev.off()
