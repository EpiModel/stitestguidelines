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
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))

haz.sti <- as.numeric(colMeans(tail(sim.base$epi$ir100.sti, 52)))
ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
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
gc.tx <- rep(NA, length(sims))
gc.txyr1 <- rep(NA, length(sims))
gc.tx.asympt <- rep(NA, length(sims))
gc.tx.asymptyr1 <- rep(NA, length(sims))
gc.tx.sympt <- rep(NA, length(sims))
gc.tx.symptyr1 <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.tx <- rep(NA, length(sims))
ct.txyr1 <- rep(NA, length(sims))
ct.tx.asympt <- rep(NA, length(sims))
ct.tx.asymptyr1 <- rep(NA, length(sims))
ct.tx.sympt <- rep(NA, length(sims))
ct.tx.symptyr1 <- rep(NA, length(sims))

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

df <- data.frame(anncov, hrcov, annint, hrint,
                 gc.incid, gc.pia, gc.txyr1, gc.tx,
                 gc.tx.asymptyr1, gc.tx.asympt, gc.tx.symptyr1, gc.tx.sympt,
                 ct.incid, ct.pia, ct.txyr1, ct.tx,
                 ct.tx.asymptyr1, ct.tx.asympt, ct.tx.symptyr1, ct.tx.sympt,
                 syph.incid, syph.pia,
                 syph.tx.earlyyr1, syph.tx.early,
                 syph.tx.lateyr1, syph.tx.late,
                 syph.txyr1, syph.tx,
                 syph.tx.asymptyr1, syph.tx.asympt,
                 syph.tx.symptyr1, syph.tx.sympt)

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

   #PIA (Cumulative)
   ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
   vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
   vec.pia.gc <- vec.nia.gc/ir.base.gc
   vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
   df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
   vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
   vec.pia.ct <- vec.nia.ct/ir.base.ct
   vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
   df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

   ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
   vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
   vec.pia.syph <- vec.nia.syph/ir.base.syph
   vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
   df$syph.pia[i] <- paste0(round(quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
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


   cat("*")

}

df

write.csv(df, "analysis/STD Table 1.csv")
