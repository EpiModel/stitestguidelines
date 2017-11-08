## STI Testing Guidelines Table 1
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
#load("data/followup/sim.n3000.rda")
load("data/sim.n3000.rda")
sim.base <- sim
#epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
ir.base.gc.g1 <- unname(colMeans(sim.base$epi$ir100.gc.tttraj1)) * 1000
ir.base.gc.g2 <- unname(colMeans(sim.base$epi$ir100.gc.tttraj2)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
ir.base.ct.g1 <- unname(colMeans(sim.base$epi$ir100.ct.tttraj1)) * 1000
ir.base.ct.g2 <- unname(colMeans(sim.base$epi$ir100.ct.tttraj2)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
ir.base.syph.g1 <- unname(colMeans(sim.base$epi$ir100.syph.tttraj1)) * 1000
ir.base.syph.g2 <- unname(colMeans(sim.base$epi$ir100.syph.tttraj2)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))

haz.sti <- as.numeric(colMeans(tail(sim.base$epi$ir100.sti, 52)))
ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
ir.base.sti.g1 <- unname(colMeans(sim.base$epi$ir100.sti.tttraj1)) * 1000
ir.base.sti.g2 <- unname(colMeans(sim.base$epi$ir100.sti.tttraj2)) * 1000
incid.base.sti <- unname(colSums(sim.base$epi$incid.sti))

# Varying Lower-Risk Coverage
# 3002, 3004, 3006, 3008
# : Annual = 10%, 20%, 30%, 40% increase, 364 days, HR = 0%, 182 days
# Varying Higher-Risk Coverage
# 3018, 3036, 3054, 3072, 3090, 3108, 3126, 3144, 3162, 3180
#  Higher-risk = 0.1 - 1.0 by 0.1, 182 days, Ann = 10%, 364 days

# Newer way:

sims <- c(3000, 3002, 3004, 3006, 3008,
          3018, 3036, 3054, 3072, 3090,
          3108, 3126, 3144, 3162, 3180)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
gc.nia <- rep(NA, length(sims))
gc.tx <- rep(NA, length(sims))
gc.txyr1 <- rep(NA, length(sims))
gc.tx.asympt <- rep(NA, length(sims))
gc.tx.asymptyr1 <- rep(NA, length(sims))
gc.tx.sympt <- rep(NA, length(sims))
gc.tx.symptyr1 <- rep(NA, length(sims))
gc.txperpy <- rep(NA, length(sims))
gc.nnt <- rep(NA, length(sims))
gc.nnt.g1 <- rep(NA, length(sims))
gc.nnt.g2 <- rep(NA, length(sims))
gc.tests <- rep(NA, length(sims))
gc.tests.py <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.nia <- rep(NA, length(sims))
ct.tx <- rep(NA, length(sims))
ct.txyr1 <- rep(NA, length(sims))
ct.tx.asympt <- rep(NA, length(sims))
ct.tx.asymptyr1 <- rep(NA, length(sims))
ct.tx.sympt <- rep(NA, length(sims))
ct.tx.symptyr1 <- rep(NA, length(sims))
ct.txperpy <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))
ct.nnt.g1 <- rep(NA, length(sims))
ct.nnt.g2 <- rep(NA, length(sims))
ct.tests <- rep(NA, length(sims))
ct.tests.py <- rep(NA, length(sims))

syph.incid <- rep(NA, length(sims))
syph.pia <- rep(NA, length(sims))
syph.nia <- rep(NA, length(sims))
syph.tx.early <- rep(NA, length(sims))
syph.tx.earlyyr1 <- rep(NA, length(sims))
syph.tx.late <- rep(NA, length(sims))
syph.tx.lateyr1 <- rep(NA, length(sims))
syph.tx <- rep(NA, length(sims))
syph.txyr1 <- rep(NA, length(sims))
syph.tx.asympt <- rep(NA, length(sims))
syph.tx.asymptyr1 <- rep(NA, length(sims))
syph.tx.sympt <- rep(NA, length(sims))
syph.tx.symptyr1 <- rep(NA, length(sims))
syph.txperpy <- rep(NA, length(sims))
syph.nnt <- rep(NA, length(sims))
syph.nnt.g1 <- rep(NA, length(sims))
syph.nnt.g2 <- rep(NA, length(sims))
syph.tests.py <- rep(NA, length(sims))
syph.tests <- rep(NA, length(sims))

sti.incid <- rep(NA, length(sims))
sti.pia <- rep(NA, length(sims))
sti.nia <- rep(NA, length(sims))
sti.tx <- rep(NA, length(sims))
sti.txyr1 <- rep(NA, length(sims))
sti.tx.asympt <- rep(NA, length(sims))
sti.tx.asymptyr1 <- rep(NA, length(sims))
sti.tx.sympt <- rep(NA, length(sims))
sti.tx.symptyr1 <- rep(NA, length(sims))
sti.txperpy <- rep(NA, length(sims))
sti.nnt <- rep(NA, length(sims))
sti.nnt.g1 <- rep(NA, length(sims))
sti.nnt.g2 <- rep(NA, length(sims))
sti.nnt.stand <- rep(NA, length(sims))
sti.tests.py <- rep(NA, length(sims))
sti.tests <- rep(NA, length(sims))

median <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint,
                 gc.incid, gc.nia, gc.pia, gc.txyr1, gc.tx, gc.nnt, gc.nnt.g1, gc.nnt.g2,
                 gc.tx.asymptyr1, gc.tx.asympt, gc.tx.symptyr1, gc.tx.sympt, gc.txperpy, gc.tests.py, gc.tests,
                 ct.incid, ct.nia, ct.pia, ct.txyr1, ct.tx, ct.nnt, ct.nnt.g1, ct.nnt.g2,
                 ct.tx.asymptyr1, ct.tx.asympt, ct.tx.symptyr1, ct.tx.sympt, ct.txperpy, ct.tests.py, ct.tests,
                 syph.incid, syph.nia, syph.pia,
                 syph.tx.earlyyr1, syph.tx.early,
                 syph.tx.lateyr1, syph.tx.late,
                 syph.txyr1, syph.tx, syph.nnt, syph.nnt.g1, syph.nnt.g2,
                 syph.tx.asymptyr1, syph.tx.asympt,
                 syph.tx.symptyr1, syph.tx.sympt, syph.txperpy, syph.tests.py,
                 sti.incid, sti.nia, sti.pia, sti.txyr1, sti.tx, sti.nnt, sti.nnt.stand, sti.nnt.g1, sti.nnt.g2,
                 sti.tx.asymptyr1, sti.tx.asympt, sti.tx.symptyr1, sti.tx.sympt,
                 sti.txperpy, sti.tests.py, sti.tests,
                 median)

for (i in seq_along(sims)) {

   fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
   load(fn)

   df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
   df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
   df$annint[i] <- sim$param$stitest.active.int
   df$hrint[i] <- sim$param$sti.highrisktest.int

   # Incidence Rate over the first year
   vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
   vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
   vec.ir.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
   vec.ir.sti <- unname(colMeans(tail(sim$epi$ir100.sti, 52)))

   df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

   df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
   df$syph.incid[i] <- paste0(round(quantile(vec.ir.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
   df$sti.incid[i] <- paste0(round(quantile(vec.ir.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

   #PIA (Cumulative)
   ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
   vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
   vec.pia.gc <- vec.nia.gc/ir.base.gc
   vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
   df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
   df$gc.nia[i] <- paste0(round(quantile(vec.nia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.nia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.nia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
   vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
   vec.pia.ct <- vec.nia.ct/ir.base.ct
   vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
   df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
   df$ct.nia[i] <- paste0(round(quantile(vec.nia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.nia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.nia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
   vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
   vec.pia.syph <- vec.nia.syph/ir.base.syph
   vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
   df$syph.pia[i] <- paste0(round(quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
   df$syph.nia[i] <- paste0(round(quantile(vec.nia.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.nia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.nia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti)) * 1000
   vec.nia.sti <- round(ir.base.sti - ir.comp.sti, 1)
   vec.pia.sti <- vec.nia.sti/ir.base.sti
   vec.pia.sti <- vec.pia.sti[vec.pia.sti > -Inf]
   df$sti.pia[i] <- paste0(round(quantile(vec.pia.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
   df$sti.nia[i] <- paste0(round(quantile(vec.nia.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.nia.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.nia.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  # Tx
   vec.txCT <- unname(colSums(sim$epi$txCT))
   vec.txCTyr1 <- as.numeric(colSums(head(sim$epi$txCT, 52)))
   vec.txCT.asympt <- unname(colSums(sim$epi$txCT_asympt))
   vec.txCT.asympt.yr1 <- as.numeric(colSums(head(sim$epi$txCT_asympt, 52)))
   vec.txCT.sympt <- unname(colSums(sim$epi$txCT - sim$epi$txCT_asympt))
   vec.txCT.sympt.yr1 <- as.numeric(colSums(head(sim$epi$txCT - sim$epi$txCT_asympt, 52)))
   df$ct.tx[i] <- paste0(round(quantile(vec.txCT, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.txCT, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.txCT, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
   df$ct.txyr1[i] <- paste0(round(quantile(vec.txCTyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.txCTyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.txCTyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
   df$ct.tx.asympt[i] <- paste0(round(quantile(vec.txCT.asympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.txCT.asympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.txCT.asympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
   df$ct.tx.asymptyr1[i] <- paste0(round(quantile(vec.txCT.asympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.txCT.asympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.txCT.asympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
   df$ct.tx.sympt[i] <- paste0(round(quantile(vec.txCT.sympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(vec.txCT.sympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(vec.txCT.sympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")
   df$ct.tx.symptyr1[i] <- paste0(round(quantile(vec.txCT.sympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                   " (", round(quantile(vec.txCT.sympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                   " - ", round(quantile(vec.txCT.sympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                   ")")


   vec.txGC <- unname(colSums(sim$epi$txGC))
   vec.txGCyr1 <- as.numeric(colSums(head(sim$epi$txGC, 52)))
   vec.txGC.asympt <- unname(colSums(sim$epi$txGC_asympt))
   vec.txGC.asympt.yr1 <- as.numeric(colSums(head(sim$epi$txGC_asympt, 52)))
   vec.txGC.sympt <- unname(colSums(sim$epi$txGC - sim$epi$txGC_asympt))
   vec.txGC.sympt.yr1 <- as.numeric(colSums(head(sim$epi$txGC - sim$epi$txGC_asympt, 52)))
   df$gc.tx[i] <- paste0(round(quantile(vec.txGC, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.txGC, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.txGC, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
   df$gc.txyr1[i] <- paste0(round(quantile(vec.txGCyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.txGCyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.txGCyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
   df$gc.tx.asympt[i] <- paste0(round(quantile(vec.txGC.asympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(vec.txGC.asympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(vec.txGC.asympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")
   df$gc.tx.asymptyr1[i] <- paste0(round(quantile(vec.txGC.asympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                   " (", round(quantile(vec.txGC.asympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                   " - ", round(quantile(vec.txGC.asympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                   ")")
   df$gc.tx.sympt[i] <- paste0(round(quantile(vec.txGC.sympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(vec.txGC.sympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(vec.txGC.sympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$gc.tx.symptyr1[i] <- paste0(round(quantile(vec.txGC.sympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                  " (", round(quantile(vec.txGC.sympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                  " - ", round(quantile(vec.txGC.sympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                  ")")

   vec.txearlysyph <- unname(colSums(sim$epi$txearlysyph))
   vec.txearlysyphyr1 <- as.numeric(colSums(head(sim$epi$txearlysyph, 52)))
   df$syph.tx.early[i] <- paste0(round(quantile(vec.txearlysyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(vec.txearlysyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(vec.txearlysyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$syph.tx.earlyyr1[i] <- paste0(round(quantile(vec.txearlysyphyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                 " (", round(quantile(vec.txearlysyphyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                 " - ", round(quantile(vec.txearlysyphyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                 ")")

   vec.txlatesyph <- unname(colSums(sim$epi$txlatesyph))
   vec.txlatesyphyr1 <- as.numeric(colSums(head(sim$epi$txlatesyph, 52)))
   df$syph.tx.late[i] <- paste0(round(quantile(vec.txlatesyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                              " (", round(quantile(vec.txlatesyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                              " - ", round(quantile(vec.txlatesyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                              ")")
   df$syph.tx.lateyr1[i] <- paste0(round(quantile(vec.txlatesyphyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(vec.txlatesyphyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(vec.txlatesyphyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")

   vec.txsyph <- unname(colSums(sim$epi$txsyph))
   vec.txsyphyr1 <- as.numeric(colSums(head(sim$epi$txearlysyph, 52)))
   df$syph.txyr1[i] <- paste0(round(quantile(vec.txsyphyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(vec.txsyphyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(vec.txsyphyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
   df$syph.tx[i] <- paste0(round(quantile(vec.txsyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                 " (", round(quantile(vec.txsyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                 " - ", round(quantile(vec.txsyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                 ")")

   vec.txsyph.asympt <- unname(colSums(sim$epi$txsyph_asympt))
   vec.txsyph.asympt.yr1 <- as.numeric(colSums(head(sim$epi$txsyph_asympt, 52)))
   vec.txsyph.sympt <- unname(colSums(sim$epi$txsyph - sim$epi$txsyph_asympt))
   vec.txsyph.sympt.yr1 <- as.numeric(colSums(head(sim$epi$txsyph - sim$epi$txsyph_asympt, 52)))

   df$syph.tx.asympt[i] <- paste0(round(quantile(vec.txsyph.asympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(vec.txsyph.asympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(vec.txsyph.asympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")
   df$syph.tx.asymptyr1[i] <- paste0(round(quantile(vec.txsyph.asympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                   " (", round(quantile(vec.txsyph.asympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                   " - ", round(quantile(vec.txsyph.asympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                   ")")
   df$syph.tx.sympt[i] <- paste0(round(quantile(vec.txsyph.sympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(vec.txsyph.sympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(vec.txsyph.sympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$syph.tx.symptyr1[i] <- paste0(round(quantile(vec.txsyph.sympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                  " (", round(quantile(vec.txsyph.sympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                  " - ", round(quantile(vec.txsyph.sympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                  ")")


   vec.txSTI <- unname(colSums(sim$epi$txSTI))
   vec.txSTIyr1 <- as.numeric(colSums(head(sim$epi$txSTI, 52)))
   vec.txSTI.asympt <- unname(colSums(sim$epi$txSTI_asympt))
   vec.txSTI.asympt.yr1 <- as.numeric(colSums(head(sim$epi$txSTI_asympt, 52)))
   vec.txSTI.sympt <- unname(colSums(sim$epi$txSTI - sim$epi$txSTI_asympt))
   vec.txSTI.sympt.yr1 <- as.numeric(colSums(head(sim$epi$txSTI - sim$epi$txSTI_asympt, 52)))
   df$sti.tx[i] <- paste0(round(quantile(vec.txSTI, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.txSTI, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.txSTI, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
   df$sti.txyr1[i] <- paste0(round(quantile(vec.txSTIyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.txSTIyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.txSTIyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
   df$sti.tx.asympt[i] <- paste0(round(quantile(vec.txSTI.asympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(vec.txSTI.asympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(vec.txSTI.asympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")
   df$sti.tx.asymptyr1[i] <- paste0(round(quantile(vec.txSTI.asympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                   " (", round(quantile(vec.txSTI.asympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                   " - ", round(quantile(vec.txSTI.asympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                   ")")
   df$sti.tx.sympt[i] <- paste0(round(quantile(vec.txSTI.sympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(vec.txSTI.sympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(vec.txSTI.sympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$sti.tx.symptyr1[i] <- paste0(round(quantile(vec.txSTI.sympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                  " (", round(quantile(vec.txSTI.sympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                  " - ", round(quantile(vec.txSTI.sympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                  ")")

   # Tests
   gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
   gc.sympt.tests <- unname(colSums(sim$epi$GCsympttests, na.rm = TRUE))
   gc.tests <- gc.asympt.tests + gc.sympt.tests

   gc.asympt.tests.g1 <- unname(colSums(sim$epi$GCasympttests.tttraj1, na.rm = TRUE))
   gc.sympt.tests.g1 <- unname(colSums(sim$epi$GCsympttests.tttraj1, na.rm = TRUE))
   gc.tests.g1 <- gc.asympt.tests.g1 + gc.sympt.tests.g1

   gc.asympt.tests.g2 <- unname(colSums(sim$epi$GCasympttests.tttraj2, na.rm = TRUE))
   gc.sympt.tests.g2 <- unname(colSums(sim$epi$GCsympttests.tttraj2, na.rm = TRUE))
   gc.tests.g2 <- gc.asympt.tests.g2 + gc.sympt.tests.g2

   ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
   ct.sympt.tests <- unname(colSums(sim$epi$CTsympttests, na.rm = TRUE))
   ct.tests <- ct.asympt.tests + ct.sympt.tests

   ct.asympt.tests.g1 <- unname(colSums(sim$epi$CTasympttests.tttraj1, na.rm = TRUE))
   ct.sympt.tests.g1 <- unname(colSums(sim$epi$CTsympttests.tttraj1, na.rm = TRUE))
   ct.tests.g1 <- ct.asympt.tests.g1 + ct.sympt.tests.g1

   ct.asympt.tests.g2 <- unname(colSums(sim$epi$CTasympttests.tttraj2, na.rm = TRUE))
   ct.sympt.tests.g2 <- unname(colSums(sim$epi$CTsympttests.tttraj2, na.rm = TRUE))
   ct.tests.g2 <- ct.asympt.tests.g2 + ct.sympt.tests.g2

   syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
   syph.sympt.tests <- unname(colSums(sim$epi$syphsympttests, na.rm = TRUE))
   syph.tests <- syph.asympt.tests + syph.sympt.tests

   syph.asympt.tests.g1 <- unname(colSums(sim$epi$syphasympttests.tttraj1, na.rm = TRUE))
   syph.sympt.tests.g1 <- unname(colSums(sim$epi$syphsympttests.tttraj1, na.rm = TRUE))
   syph.tests.g1 <- syph.asympt.tests.g1 + syph.sympt.tests.g1

   syph.asympt.tests.g2 <- unname(colSums(sim$epi$syphasympttests.tttraj2, na.rm = TRUE))
   syph.sympt.tests.g2 <- unname(colSums(sim$epi$syphsympttests.tttraj2, na.rm = TRUE))
   syph.tests.g2 <- syph.asympt.tests.g2 + syph.sympt.tests.g2

   sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))
   sti.sympt.tests <- unname(colSums(sim$epi$stisympttests, na.rm = TRUE))
   sti.tests <- sti.asympt.tests + sti.sympt.tests

   py <- unname(colSums(sim$epi$num, na.rm = TRUE))

   gc.tests.py <- 52 * (gc.tests / py)
   ct.tests.py <-  52 * (ct.tests / py)
   syph.tests.py <-  52 * (syph.tests / py)
   sti.tests.py <-  52 * (sti.tests / py)

   df$gc.tests[i] <- paste0(round(quantile(gc.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(gc.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(gc.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$ct.tests[i] <- paste0(round(quantile(ct.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(ct.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(ct.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$syph.tests[i] <- paste0(round(quantile(syph.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                 " (", round(quantile(syph.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                 " - ", round(quantile(syph.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                 ")")
   df$sti.tests[i] <- paste0(round(quantile(sti.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(sti.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(sti.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")

   df$gc.tests.py[i] <- paste0(round(quantile(gc.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(gc.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(gc.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
   df$ct.tests.py[i] <- paste0(round(quantile(ct.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(ct.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(ct.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
   df$syph.tests.py[i] <- paste0(round(quantile(syph.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(syph.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(syph.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
   df$sti.tests.py[i] <- paste0(round(quantile(sti.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(sti.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(sti.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")

   # Tx per person-year infected
   vec.tx.gcpy <- unname(colMeans(52 * sim$epi$txGC / (sim$epi$num * sim$epi$prev.gc)))
   df$gc.txperpy[i] <- paste0(round(quantile(vec.tx.gcpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.tx.gcpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.tx.gcpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
   vec.tx.ctpy <- unname(colMeans(52 * sim$epi$txCT / (sim$epi$num * sim$epi$prev.ct)))
   df$ct.txperpy[i] <- paste0(round(quantile(vec.tx.ctpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.tx.ctpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.tx.ctpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
   vec.tx.syphpy <- unname(colMeans(52 * sim$epi$txsyph / (sim$epi$num * sim$epi$prev.syph)))
   df$syph.txperpy[i] <- paste0(round(quantile(vec.tx.syphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.tx.syphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.tx.syphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
   vec.tx.stipy <- unname(colMeans(52 * sim$epi$txSTI / (sim$epi$num * sim$epi$prev.sti)))
   df$sti.txperpy[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

   # Number needed to treat
   gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
   gc.asympt.tests.g1 <- unname(colSums(sim$epi$GCasympttests.tttraj1, na.rm = TRUE))
   gc.asympt.tests.g2 <- unname(colSums(sim$epi$GCasympttests.tttraj2, na.rm = TRUE))
   ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
   ct.asympt.tests.g1 <- unname(colSums(sim$epi$CTasympttests.tttraj1, na.rm = TRUE))
   ct.asympt.tests.g2 <- unname(colSums(sim$epi$CTasympttests.tttraj2, na.rm = TRUE))
   syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
   syph.asympt.tests.g1 <- unname(colSums(sim$epi$syphasympttests.tttraj1, na.rm = TRUE))
   syph.asympt.tests.g2 <- unname(colSums(sim$epi$syphasympttests.tttraj2, na.rm = TRUE))
   sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))
   sti.asympt.tests.g1 <- unname(colSums(sim$epi$stiasympttests.tttraj1, na.rm = TRUE))
   sti.asympt.tests.g2 <- unname(colSums(sim$epi$stiasympttests.tttraj2, na.rm = TRUE))

   vec.gc.nnt <- gc.asympt.tests / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
   vec.gc.nnt.g1 <- gc.asympt.tests.g1 / (median(ir.base.gc.g1) - (unname(colSums(sim$epi$ir100.gc.tttraj1)) * 1000))
   vec.gc.nnt.g2 <- gc.asympt.tests.g2 / (median(ir.base.gc.g2) - (unname(colSums(sim$epi$ir100.gc.tttraj2)) * 1000))
   vec.ct.nnt <- ct.asympt.tests / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))
   vec.ct.nnt.g1 <- ct.asympt.tests.g1 / (median(ir.base.ct.g1) - (unname(colSums(sim$epi$ir100.ct.tttraj1)) * 1000))
   vec.ct.nnt.g2 <- ct.asympt.tests.g2 / (median(ir.base.ct.g2) - (unname(colSums(sim$epi$ir100.ct.tttraj2)) * 1000))
   vec.syph.nnt <- syph.asympt.tests / (median(incid.base.syph) - unname(colSums(sim$epi$incid.syph)))
   vec.syph.nnt.g1 <- syph.asympt.tests.g1 / (median(ir.base.syph.g1) - (unname(colSums(sim$epi$ir100.syph.tttraj1)) * 1000))
   vec.syph.nnt.g2 <- syph.asympt.tests.g2 / (median(ir.base.syph.g2) - (unname(colSums(sim$epi$ir100.syph.tttraj2)) * 1000))
   vec.sti.nnt <- sti.asympt.tests / (median(incid.base.sti) - unname(colSums(sim$epi$incid.sti)))
   vec.sti.nnt.g1 <- sti.asympt.tests.g1 / (median(ir.base.sti.g1) - (unname(colSums(sim$epi$ir100.sti.tttraj1)) * 1000))
   vec.sti.nnt.g2 <- sti.asympt.tests.g2 / (median(ir.base.sti.g2) - (unname(colSums(sim$epi$ir100.sti.tttraj2)) * 1000))

   vec.sti.nnt.stand <- sti.asympt.tests / (median(vec.nia.sti))

   df$sti.nnt.stand[i] <- paste0(round(quantile(vec.sti.nnt.stand, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                 " (", round(quantile(vec.sti.nnt.stand, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                 " - ", round(quantile(vec.sti.nnt.stand, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                 ")")


   df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
   df$gc.nnt.g1[i] <- paste0(round(quantile(vec.gc.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.gc.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.gc.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
   df$gc.nnt.g2[i] <- paste0(round(quantile(vec.gc.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.gc.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.gc.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
   df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
   df$ct.nnt.g1[i] <- paste0(round(quantile(vec.ct.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ct.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ct.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
   df$ct.nnt.g2[i] <- paste0(round(quantile(vec.ct.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ct.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ct.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
   df$syph.nnt[i] <- paste0(round(quantile(vec.syph.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.syph.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.syph.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
   df$syph.nnt.g1[i] <- paste0(round(quantile(vec.syph.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.syph.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.syph.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
   df$syph.nnt.g2[i] <- paste0(round(quantile(vec.syph.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.syph.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.syph.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
   df$sti.nnt[i] <- paste0(round(quantile(vec.sti.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.sti.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.sti.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
   df$sti.nnt.g1[i] <- paste0(round(quantile(vec.sti.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.sti.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.sti.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
   df$sti.nnt.g2[i] <- paste0(round(quantile(vec.sti.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.sti.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.sti.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

   meddiff <- (median(incid.base.sti) - unname(colSums(sim$epi$incid.sti)))
   df$median[i] <- paste0(round(quantile(meddiff, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(meddiff, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(meddiff, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

    cat("*")

}

df

write.csv(df, "analysis/STD Table 1.csv")
