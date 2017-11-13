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

ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))
tests.base <- unname(colSums(sim.base$epi$hivtests.nprep))

ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
ir.base.gc.g1 <- unname(colMeans(sim.base$epi$ir100.gc.tttraj1)) * 1000
ir.base.gc.g2 <- unname(colMeans(sim.base$epi$ir100.gc.tttraj2)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))

ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
ir.base.ct.g1 <- unname(colMeans(sim.base$epi$ir100.ct.tttraj1)) * 1000
ir.base.ct.g2 <- unname(colMeans(sim.base$epi$ir100.ct.tttraj2)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))

ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
ir.base.syph.g1 <- unname(colMeans(sim.base$epi$ir100.syph.tttraj1)) * 1000
ir.base.syph.g2 <- unname(colMeans(sim.base$epi$ir100.syph.tttraj2)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests))

ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
ir.base.sti.g1 <- unname(colMeans(sim.base$epi$ir100.sti.tttraj1)) * 1000
ir.base.sti.g2 <- unname(colMeans(sim.base$epi$ir100.sti.tttraj2)) * 1000
incid.base.sti <- unname(colSums(sim.base$epi$incid.sti))
tests.sti.base <- unname(colSums(sim.base$epi$stiasympttests))

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
gc.tx <- rep(NA, length(sims))
gc.txyr1 <- rep(NA, length(sims))
gc.tx.asympt <- rep(NA, length(sims))
gc.tx.asymptyr1 <- rep(NA, length(sims))
gc.tx.sympt <- rep(NA, length(sims))
gc.tx.symptyr1 <- rep(NA, length(sims))
gc.txperpy <- rep(NA, length(sims))
gc.nnt <- rep(NA, length(sims))
gc.asympt.tests <- rep(NA, length(sims))
gc.asympt.tests.py <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.tx <- rep(NA, length(sims))
ct.txyr1 <- rep(NA, length(sims))
ct.tx.asympt <- rep(NA, length(sims))
ct.tx.asymptyr1 <- rep(NA, length(sims))
ct.tx.sympt <- rep(NA, length(sims))
ct.tx.symptyr1 <- rep(NA, length(sims))
ct.txperpy <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))
ct.asympt.tests <- rep(NA, length(sims))
ct.asympt.tests.py <- rep(NA, length(sims))

syph.incid <- rep(NA, length(sims))
syph.pia <- rep(NA, length(sims))
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
syph.earlytxperpy <- rep(NA, length(sims))
syph.latetxperpy <- rep(NA, length(sims))
syph.nnt <- rep(NA, length(sims))
syph.asympt.tests.py <- rep(NA, length(sims))
syph.asympt.tests <- rep(NA, length(sims))

sti.incid <- rep(NA, length(sims))
sti.pia <- rep(NA, length(sims))
sti.tx <- rep(NA, length(sims))
sti.txyr1 <- rep(NA, length(sims))
sti.tx.asympt <- rep(NA, length(sims))
sti.tx.asymptyr1 <- rep(NA, length(sims))
sti.tx.sympt <- rep(NA, length(sims))
sti.tx.symptyr1 <- rep(NA, length(sims))
sti.txperpy <- rep(NA, length(sims))
sti.nnt <- rep(NA, length(sims))
sti.asympt.tests.py <- rep(NA, length(sims))
sti.asympt.tests <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint,
                 gc.incid, gc.pia, gc.txyr1, gc.tx, gc.nnt,
                 gc.tx.asymptyr1, gc.tx.asympt, gc.tx.symptyr1, gc.tx.sympt,
                 gc.txperpy, gc.asympt.tests.py, gc.asympt.tests,
                 ct.incid, ct.pia, ct.txyr1, ct.tx, ct.nnt,
                 ct.tx.asymptyr1, ct.tx.asympt, ct.tx.symptyr1, ct.tx.sympt,
                 ct.txperpy, ct.asympt.tests.py, ct.asympt.tests,
                 syph.incid, syph.pia,
                 syph.tx.earlyyr1, syph.tx.early,
                 syph.tx.lateyr1, syph.tx.late,
                 syph.txyr1, syph.tx, syph.nnt,
                 syph.tx.asymptyr1, syph.tx.asympt,
                 syph.tx.symptyr1, syph.tx.sympt, syph.txperpy,
                 syph.earlytxperpy, syph.latetxperpy, syph.asympt.tests.py, syph.asympt.tests,
                 sti.incid, sti.pia, sti.txyr1, sti.tx, sti.nnt,
                 sti.tx.asymptyr1, sti.tx.asympt, sti.tx.symptyr1, sti.tx.sympt,
                 sti.txperpy, sti.asympt.tests.py, sti.asympt.tests)

for (i in seq_along(sims)) {

   #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
   fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
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
   incid.gc <- unname(colSums(sim$epi$incid.gc))
   vec.nia.gc <- incid.base.gc - incid.gc
   vec.pia.gc <- vec.nia.gc/incid.base.gc
   df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   incid.ct <- unname(colSums(sim$epi$incid.ct))
   vec.nia.ct <- incid.base.ct - incid.ct
   vec.pia.ct <- vec.nia.ct/incid.base.ct
   df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   incid.syph <- unname(colSums(sim$epi$incid.syph))
   vec.nia.syph <- incid.base.syph - incid.syph
   vec.pia.syph <- vec.nia.syph/incid.base.syph
   df$syph.pia[i] <- paste0(round(quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

   incid.sti <- unname(colSums(sim$epi$incid.sti))
   vec.nia.sti <- incid.base.sti - incid.sti
   vec.pia.sti <- vec.nia.sti/incid.base.sti
   df$sti.pia[i] <- paste0(round(quantile(vec.pia.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
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

   gc.asympt.tests.py <- 52 * (gc.asympt.tests / py)
   ct.asympt.tests.py <-  52 * (ct.asympt.tests / py)
   syph.asympt.tests.py <-  52 * (syph.asympt.tests / py)
   sti.asympt.tests.py <-  52 * (sti.asympt.tests / py)

   df$gc.asympt.tests[i] <- paste0(round(quantile(gc.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(gc.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(gc.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$ct.asympt.tests[i] <- paste0(round(quantile(ct.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(ct.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(ct.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
   df$syph.asympt.tests[i] <- paste0(round(quantile(syph.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                 " (", round(quantile(syph.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                 " - ", round(quantile(syph.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                 ")")
   df$sti.asympt.tests[i] <- paste0(round(quantile(sti.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                " (", round(quantile(sti.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                " - ", round(quantile(sti.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                ")")

   df$gc.asympt.tests.py[i] <- paste0(round(quantile(gc.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(gc.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(gc.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
   df$ct.asympt.tests.py[i] <- paste0(round(quantile(ct.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(ct.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(ct.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
   df$syph.asympt.tests.py[i] <- paste0(round(quantile(syph.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(syph.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(syph.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
   df$sti.asympt.tests.py[i] <- paste0(round(quantile(sti.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(sti.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(sti.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
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


    # Prev.primsecosyph (Incub, prim, syph)
   df.prev.earlysyph <- sim$epi$txearlysyph[1:521,] / (sim$epi$num[1:521,] * sim$epi$prev.primsecosyph[1:521,])

   # Remove NaNs
   for (j in 1:ncol(df.prev.earlysyph)) {

     df.prev.earlysyph[which(is.nan(df.prev.earlysyph[, j])), j] <- 0.0

   }
   vec.tx.earlysyphpy <- unname(colMeans(52 * df.prev.earlysyph))

   #vec.tx.earlysyphpy <- unname(colMeans(52 * sim$epi$txearlysyph / (sim$epi$num * sim$epi$prev.primsecosyph)))
   df$syph.earlytxperpy[i] <- paste0(round(quantile(vec.tx.earlysyphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                    " (", round(quantile(vec.tx.earlysyphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                    " - ", round(quantile(vec.tx.earlysyphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                    ")")

   # Late syphilis prevalence = Prev.syph - Prev.primsecosyph
   df.prev.latesyph <- sim$epi$txlatesyph[1:521,] / (sim$epi$num[1:521,] * (sim$epi$prev.syph[1:521, ] - sim$epi$prev.primsecosyph[1:521,]))

   # Remove NaNs
   for (j in 1:ncol(df.prev.latesyph)) {

     df.prev.latesyph[which(is.nan(df.prev.latesyph[, j])), j] <- 0.0

   }
   vec.tx.latesyphpy <- unname(colMeans(52 * df.prev.latesyph))
   #vec.tx.latesyphpy <- unname(colMeans(52 * sim$epi$txlatesyph / (sim$epi$num * (sim$epi$prev.syph - sim$epi$prev.primsecosyph))))
   df$syph.latetxperpy[i] <- paste0(round(quantile(vec.tx.latesyphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                    " (", round(quantile(vec.tx.latesyphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                    " - ", round(quantile(vec.tx.latesyphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                    ")")

   vec.tx.stipy <- unname(colMeans(52 * sim$epi$txSTI / (sim$epi$num * sim$epi$prev.sti)))
   df$sti.txperpy[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

   # Number needed to treat
   # incremental tests / number of infections averted
   vec.gc.nnt <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - unname(colSums(sim$epi$incid.gc)))
   vec.ct.nnt <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - unname(colSums(sim$epi$incid.ct)))
   vec.syph.nnt <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - unname(colSums(sim$epi$incid.syph)))
   vec.sti.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base.sti - unname(colSums(sim$epi$incid.sti)))

   df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   df$syph.nnt[i] <- paste0(round(quantile(vec.syph.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.syph.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.syph.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

   df$sti.nnt[i] <- paste0(round(quantile(vec.sti.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.sti.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.sti.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
    cat("*")

}

df

write.csv(df, "analysis/STD Table 1.csv")
