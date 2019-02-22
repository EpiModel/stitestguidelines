## STI Testing Guidelines Supplementary Tables 13 and 14
# Treatment Adherence

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

# Base - No annual or high-risk
load("data/followup/Guidelines Paper/sim.n9000.rda")
sim.base <- sim

incid.base <- unname(colSums(sim.base$epi$incid))
tests.base <- unname(colSums(sim.base$epi$hivtests.nprep))

incid.base.gc <- unname(colSums(sim.base$epi$incid.rgc)) + unname(colSums(sim.base$epi$incid.ugc))
incid.base.gc.g1 <- unname(colSums(sim.base$epi$incid.rgc.tttraj2)) + unname(colSums(sim.base$epi$incid.ugc.tttraj2))
incid.base.gc.g2 <- unname(colSums(sim.base$epi$incid.rgc.tttraj2)) + unname(colSums(sim.base$epi$incid.ugc.tttraj2))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))
tests.gc.base.g1 <- unname(colSums(sim.base$epi$GCasympttests.tttraj1))
tests.gc.base.g2 <- unname(colSums(sim.base$epi$GCasympttests.tttraj2))

incid.base.ct <- unname(colSums(sim.base$epi$incid.rct)) + unname(colSums(sim.base$epi$incid.uct))
incid.base.ct.g1 <- unname(colSums(sim.base$epi$incid.rct.tttraj1)) + unname(colSums(sim.base$epi$incid.uct.tttraj1))
incid.base.ct.g2 <- unname(colSums(sim.base$epi$incid.rct.tttraj2)) + unname(colSums(sim.base$epi$incid.uct.tttraj2))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))
tests.ct.base.g1 <- unname(colSums(sim.base$epi$CTasympttests.tttraj1))
tests.ct.base.g2 <- unname(colSums(sim.base$epi$CTasympttests.tttraj2))

incid.base.gcct <- incid.base.ct + incid.base.gc
incid.base.gcct.g1 <- incid.base.gc.g1 + incid.base.ct.g1
incid.base.gcct.g2 <- incid.base.gc.g2 + incid.base.ct.g2
tests.sti.base <- tests.gc.base + tests.ct.base
tests.sti.base.g1 <- tests.gc.base.g1 + tests.ct.base.g1
tests.sti.base.g2 <- tests.ct.base.g1 + tests.ct.base.g2

# Treatment Adherence
# Baseline: 9000 (No HR - 100%)
# 100% to 0%: 9416, 9415, 9414, 9413, 9412, 9411, 9410, 9409, 9408, 9407 (5% HR)

# Newer way:
sims <- c(9000, 9407:9416)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
txprob <- rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
gc.tx <- rep(NA, length(sims))
gc.txyr1 <- rep(NA, length(sims))
gc.tx.asympt <- rep(NA, length(sims))
gc.tx.asymptyr1 <- rep(NA, length(sims))
gc.tx.sympt <- rep(NA, length(sims))
gc.tx.symptyr1 <- rep(NA, length(sims))
gctxpy <- rep(NA, length(sims))
gc.nnt <- rep(NA, length(sims))
gc.asympt.tests <- rep(NA, length(sims))
gc.asympt.tests.py <- rep(NA, length(sims))

gc.incid.g1 <- rep(NA, length(sims))
gc.pia.g1 <- rep(NA, length(sims))
gc.asympt.tests.py.g1 <- rep(NA, length(sims))
gc.asympt.tests.g1 <- rep(NA, length(sims))
gctx.g1 <- rep(NA, length(sims))
gctxpy.g1 <- rep(NA, length(sims))

gc.incid.g2 <- rep(NA, length(sims))
gc.pia.g2 <- rep(NA, length(sims))
gc.asympt.tests.py.g2 <- rep(NA, length(sims))
gc.asympt.tests.g2 <- rep(NA, length(sims))
gctx.g2 <- rep(NA, length(sims))
gctxpy.g2 <- rep(NA, length(sims))

gc.nnt.g1 <- rep(NA, length(sims))
gc.nnt.g2 <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.tx <- rep(NA, length(sims))
ct.txyr1 <- rep(NA, length(sims))
ct.tx.asympt <- rep(NA, length(sims))
ct.tx.asymptyr1 <- rep(NA, length(sims))
ct.tx.sympt <- rep(NA, length(sims))
ct.tx.symptyr1 <- rep(NA, length(sims))
cttxpy <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))
ct.asympt.tests <- rep(NA, length(sims))
ct.asympt.tests.py <- rep(NA, length(sims))

ct.incid.g1 <- rep(NA, length(sims))
ct.pia.g1 <- rep(NA, length(sims))
ct.asympt.tests.py.g1 <- rep(NA, length(sims))
ct.asympt.tests.g1 <- rep(NA, length(sims))
cttx.g1 <- rep(NA, length(sims))
cttxpy.g1 <- rep(NA, length(sims))

ct.incid.g2 <- rep(NA, length(sims))
ct.pia.g2 <- rep(NA, length(sims))
ct.asympt.tests.py.g2 <- rep(NA, length(sims))
ct.asympt.tests.g2 <- rep(NA, length(sims))
cttx.g2 <- rep(NA, length(sims))
cttxpy.g2 <- rep(NA, length(sims))

ct.nnt.g1 <- rep(NA, length(sims))
ct.nnt.g2 <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
sti.tx <- rep(NA, length(sims))
sti.txyr1 <- rep(NA, length(sims))
sti.tx.asympt <- rep(NA, length(sims))
sti.tx.asymptyr1 <- rep(NA, length(sims))
sti.tx.sympt <- rep(NA, length(sims))
sti.tx.symptyr1 <- rep(NA, length(sims))
txperpy.sti <- rep(NA, length(sims))
sti.nnt <- rep(NA, length(sims))
sti.asympt.tests.py <- rep(NA, length(sims))
sti.asympt.tests <- rep(NA, length(sims))

gcct.incid.g1 <- rep(NA, length(sims))
gcct.pia.g1 <- rep(NA, length(sims))
sti.asympt.tests.g1 <- rep(NA, length(sims))
sti.asympt.tests.py.g1 <- rep(NA, length(sims))
tx.sti.g1 <- rep(NA, length(sims))
txperpy.sti.g1 <- rep(NA, length(sims))

gcct.incid.g2 <- rep(NA, length(sims))
gcct.pia.g2 <- rep(NA, length(sims))
sti.asympt.tests.g2 <- rep(NA, length(sims))
sti.asympt.tests.py.g2 <- rep(NA, length(sims))
tx.sti.g2 <- rep(NA, length(sims))
txperpy.sti.g2 <- rep(NA, length(sims))

gcct.nnt <- rep(NA, length(sims))
gcct.nnt.g1 <- rep(NA, length(sims))
gcct.nnt.g2 <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, txprob,
                 gc.incid, gc.pia, gc.nnt, gctxpy, gc.asympt.tests.py,
                 gc.txyr1, gc.tx,
                 gc.tx.asymptyr1, gc.tx.asympt, gc.tx.symptyr1, gc.tx.sympt,
                 gc.asympt.tests,
                 ct.incid, ct.pia, ct.nnt, cttxpy, ct.asympt.tests.py,
                 ct.txyr1, ct.tx,
                 ct.tx.asymptyr1, ct.tx.asympt, ct.tx.symptyr1, ct.tx.sympt,
                  ct.asympt.tests,
                 gcct.incid, gcct.pia, sti.txyr1, sti.tx, gcct.nnt,
                 sti.tx.asymptyr1, sti.tx.asympt, sti.tx.symptyr1, sti.tx.sympt,
                 txperpy.sti, sti.asympt.tests.py, sti.asympt.tests,

                 # Group 1
                 gc.incid.g1, gc.pia.g1, gc.nnt.g1, gctxpy.g1,
                 gc.asympt.tests.py.g1, gc.asympt.tests.g1, gctx.g1,
                 ct.incid.g1, ct.pia.g1, ct.nnt.g1, cttxpy.g1,
                 ct.asympt.tests.py.g1, ct.asympt.tests.g1, cttx.g1,
                 gcct.incid.g1, gcct.pia.g1, gcct.nnt.g1, txperpy.sti.g1,
                 sti.asympt.tests.py.g1, sti.asympt.tests.g1, tx.sti.g1,

                 # Group 2
                 gc.incid.g2, gc.pia.g2, gc.nnt.g2, gctxpy.g2,
                 gc.asympt.tests.py.g2, gc.asympt.tests.g2, gctx.g2,
                 ct.incid.g2, ct.pia.g2, ct.nnt.g2, cttxpy.g2,
                 ct.asympt.tests.py.g2, ct.asympt.tests.g2, cttx.g2,
                 gcct.incid.g2, gcct.pia.g2, gcct.nnt.g2, txperpy.sti.g2,
                 sti.asympt.tests.py.g2, sti.asympt.tests.g2, tx.sti.g2)

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int
  df$txprob[i] <- sim$param$gc.asympt.prob.tx

  # Incidence Rate over the first year
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  vec.ir.gc.g1 <- unname(colMeans(tail(sim$epi$ir100.gc.tttraj1, 52)))
  vec.ir.gc.g2 <- unname(colMeans(tail(sim$epi$ir100.gc.tttraj2, 52)))
  df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$gc.incid.g1[i] <- paste0(round(quantile(vec.ir.gc.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.gc.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.gc.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$gc.incid.g2[i] <- paste0(round(quantile(vec.ir.gc.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.gc.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.gc.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
  vec.ir.ct.g1 <- unname(colMeans(tail(sim$epi$ir100.ct.tttraj1, 52)))
  vec.ir.ct.g2 <- unname(colMeans(tail(sim$epi$ir100.ct.tttraj2, 52)))
  df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$ct.incid.g1[i] <- paste0(round(quantile(vec.ir.ct.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.ct.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.ct.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$ct.incid.g2[i] <- paste0(round(quantile(vec.ir.ct.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.ct.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.ct.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  vec.ir.gcct <- vec.ir.gc + vec.ir.ct
  vec.ir.gcct.g1 <- vec.ir.gc.g1 + vec.ir.ct.g1
  vec.ir.gcct.g2 <- vec.ir.gc.g2 + vec.ir.ct.g2
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
  df$gcct.incid.g1[i] <- paste0(round(quantile(vec.ir.gcct.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.ir.gcct.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.ir.gcct.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
  df$gcct.incid.g2[i] <- paste0(round(quantile(vec.ir.gcct.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.ir.gcct.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.ir.gcct.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")

  #PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.rgc)) + unname(colSums(sim$epi$incid.ugc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.gc.g1 <- unname(colSums(sim$epi$incid.rgc.tttraj1)) + unname(colSums(sim$epi$incid.ugc.tttraj1))
  vec.nia.gc.g1 <- round(incid.base.gc.g1 - incid.gc.g1, 1)
  vec.pia.gc.g1 <- vec.nia.gc.g1/incid.gc.g1

  incid.gc.g2 <- unname(colSums(sim$epi$incid.rgc.tttraj2)) + unname(colSums(sim$epi$incid.ugc.tttraj2))
  vec.nia.gc.g2 <- round(incid.base.gc.g2 - incid.gc.g2, 1)
  vec.pia.gc.g2 <- vec.nia.gc.g2/incid.gc.g2

  incid.ct <- unname(colSums(sim$epi$incid.rct)) + unname(colSums(sim$epi$incid.uct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.ct.g1 <- unname(colSums(sim$epi$incid.rct.tttraj1)) + unname(colSums(sim$epi$incid.uct.tttraj1))
  vec.nia.ct.g1 <- round(incid.base.ct.g1 - incid.ct.g1, 1)
  vec.pia.ct.g1 <- vec.nia.ct.g1/incid.ct.g1

  incid.ct.g2 <- unname(colSums(sim$epi$incid.rct.tttraj2)) + unname(colSums(sim$epi$incid.uct.tttraj2))
  vec.nia.ct.g2 <- round(incid.base.ct.g2 - incid.ct.g2, 1)
  vec.pia.ct.g2 <- vec.nia.ct.g2/incid.ct.g2

  incid.gcct <- incid.gc + incid.ct
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  incid.gcct.g1 <- incid.gc.g1 + incid.ct.g1
  vec.nia.gcct.g1 <- round(incid.base.gcct.g1 - incid.gcct.g1, 1)
  vec.pia.gcct.g1 <- vec.nia.gcct.g1/incid.gcct.g1

  incid.gcct.g2 <- incid.gc.g2 + incid.ct.g2
  vec.nia.gcct.g2 <- round(incid.base.gcct.g2 - incid.gcct.g2, 1)
  vec.pia.gcct.g2 <- vec.nia.gcct.g2/incid.gcct.g2

  df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$gc.pia.g1[i] <- paste0(round(quantile(vec.pia.gc.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.gc.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.gc.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
  df$gc.pia.g2[i] <- paste0(round(quantile(vec.pia.gc.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.gc.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.gc.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
  df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$ct.pia.g1[i] <- paste0(round(quantile(vec.pia.ct.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.ct.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.ct.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
  df$ct.pia.g2[i] <- paste0(round(quantile(vec.pia.ct.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.pia.ct.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.pia.ct.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$gcct.pia.g1[i] <- paste0(round(quantile(vec.pia.gcct.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.pia.gcct.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.pia.gcct.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$gcct.pia.g2[i] <- paste0(round(quantile(vec.pia.gcct.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.pia.gcct.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.pia.gcct.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
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

  vec.txSTI <- vec.txGC + vec.txCT
  vec.txSTIyr1 <- vec.txGCyr1 + vec.txCTyr1
  vec.txSTI.asympt <- vec.txGC.asympt + vec.txCT.asympt
  vec.txSTI.asympt.yr1 <- vec.txGC.asympt.yr1 + vec.txCT.asympt.yr1
  vec.txSTI.sympt <- vec.txGC.sympt + vec.txCT.sympt
  vec.txSTI.sympt.yr1 <- vec.txGC.sympt.yr1 + vec.txCT.sympt.yr1
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

  sti.asympt.tests <- gc.asympt.tests + ct.asympt.tests
  sti.sympt.tests <- gc.sympt.tests + ct.sympt.tests
  sti.tests <- sti.asympt.tests + sti.sympt.tests

  sti.asympt.tests.g1 <- gc.asympt.tests.g1 + ct.asympt.tests.g1
  sti.sympt.tests.g1 <- gc.sympt.tests.g1 + ct.sympt.tests.g1
  sti.tests.g1 <- sti.asympt.tests.g1 + sti.sympt.tests.g1

  sti.asympt.tests.g2 <- gc.asympt.tests.g2 + ct.asympt.tests.g2
  sti.sympt.tests.g2 <- gc.sympt.tests.g2 + ct.sympt.tests.g2
  sti.tests.g2 <- sti.asympt.tests.g2 + sti.sympt.tests.g2

  py <- unname(colSums(sim$epi$num, na.rm = TRUE))
  py.g1 <- unname(colSums(sim$epi$tt.traj.sti1, na.rm = TRUE))
  py.g2 <- unname(colSums(sim$epi$tt.traj.sti2, na.rm = TRUE))

  gc.asympt.tests.py <-  52 * (gc.asympt.tests / py)
  gc.asympt.tests.py.g1 <-  52 * (gc.asympt.tests.g1 / py.g1)
  gc.asympt.tests.py.g2 <-  52 * (gc.asympt.tests.g2 / py.g2)
  ct.asympt.tests.py <-  52 * (ct.asympt.tests / py)
  ct.asympt.tests.py.g1 <-  52 * (ct.asympt.tests.g1 / py.g1)
  ct.asympt.tests.py.g2 <-  52 * (ct.asympt.tests.g2 / py.g2)
  sti.asympt.tests.py <-  52 * (sti.asympt.tests / py)
  sti.asympt.tests.py.g1 <-  52 * (sti.asympt.tests.g1 / py)
  sti.asympt.tests.py.g2 <-  52 * (sti.asympt.tests.g2 / py)

  df$gc.asympt.tests[i] <- paste0(round(quantile(gc.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                  " (", round(quantile(gc.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                  " - ", round(quantile(gc.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                  ")")
  df$gc.asympt.tests.g1[i] <- paste0(round(quantile(gc.asympt.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                     " (", round(quantile(gc.asympt.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                     " - ", round(quantile(gc.asympt.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                     ")")
  df$gc.asympt.tests.g2[i] <- paste0(round(quantile(gc.asympt.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                     " (", round(quantile(gc.asympt.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                     " - ", round(quantile(gc.asympt.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                     ")")
  df$gc.asympt.tests.py[i] <- paste0(round(quantile(gc.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(gc.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(gc.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")

  df$gc.asympt.tests.py.g1[i] <- paste0(round(quantile(gc.asympt.tests.py.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(gc.asympt.tests.py.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(gc.asympt.tests.py.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")
  df$gc.asympt.tests.py.g2[i] <- paste0(round(quantile(gc.asympt.tests.py.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(gc.asympt.tests.py.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(gc.asympt.tests.py.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")

  df$ct.asympt.tests[i] <- paste0(round(quantile(ct.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                  " (", round(quantile(ct.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                  " - ", round(quantile(ct.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                  ")")
  df$ct.asympt.tests.g1[i] <- paste0(round(quantile(ct.asympt.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                     " (", round(quantile(ct.asympt.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                     " - ", round(quantile(ct.asympt.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                     ")")
  df$ct.asympt.tests.g2[i] <- paste0(round(quantile(ct.asympt.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                     " (", round(quantile(ct.asympt.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                     " - ", round(quantile(ct.asympt.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                     ")")
  df$ct.asympt.tests.py[i] <- paste0(round(quantile(ct.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(ct.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(ct.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")

  df$ct.asympt.tests.py.g1[i] <- paste0(round(quantile(ct.asympt.tests.py.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(ct.asympt.tests.py.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(ct.asympt.tests.py.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")
  df$ct.asympt.tests.py.g2[i] <- paste0(round(quantile(ct.asympt.tests.py.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(ct.asympt.tests.py.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(ct.asympt.tests.py.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")

  df$sti.asympt.tests[i] <- paste0(round(quantile(sti.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                   " (", round(quantile(sti.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                   " - ", round(quantile(sti.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                   ")")
  df$sti.asympt.tests.g1[i] <- paste0(round(quantile(sti.asympt.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                      " (", round(quantile(sti.asympt.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                      " - ", round(quantile(sti.asympt.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                      ")")
  df$sti.asympt.tests.g2[i] <- paste0(round(quantile(sti.asympt.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                      " (", round(quantile(sti.asympt.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                      " - ", round(quantile(sti.asympt.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                      ")")
  df$sti.asympt.tests.py[i] <- paste0(round(quantile(sti.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                      " (", round(quantile(sti.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                      " - ", round(quantile(sti.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                      ")")

  df$sti.asympt.tests.py.g1[i] <- paste0(round(quantile(sti.asympt.tests.py.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                         " (", round(quantile(sti.asympt.tests.py.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                         " - ", round(quantile(sti.asympt.tests.py.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                         ")")
  df$sti.asympt.tests.py.g2[i] <- paste0(round(quantile(sti.asympt.tests.py.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                         " (", round(quantile(sti.asympt.tests.py.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                         " - ", round(quantile(sti.asympt.tests.py.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                         ")")

  # Tx per person-year infected
  vec.tx.gcpy <- unname(colMeans(52 * sim$epi$txGC / (sim$epi$num * sim$epi$prev.gc), na.rm = TRUE))
  df$gctxpy[i] <- paste0(round(quantile(vec.tx.gcpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.tx.gcpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.tx.gcpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  vec.tx.gcpy.g1 <- unname(colMeans(52 * sim$epi$txGC.tttraj1 / (sim$epi$tt.traj.gc1 * sim$epi$prev.gc.tttraj1)))
  df.prev.gc.tttraj2 <- sim$epi$txGC.tttraj2[1:521,] / (sim$epi$tt.traj.gc2[1:521,] * sim$epi$prev.gc.tttraj2[1:521,])

  # Remove NaNs
  for (j in 1:ncol(df.prev.gc.tttraj2)) {
    df.prev.gc.tttraj2[which(is.nan(df.prev.gc.tttraj2[, j])), j] <- 0.0
  }
  vec.tx.gcpy.g2 <- unname(colMeans(52 * df.prev.gc.tttraj2))
  df$gctxpy.g1[i] <- paste0(round(quantile(vec.tx.gcpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.gcpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.gcpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$gctxpy.g2[i] <- paste0(round(quantile(vec.tx.gcpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.gcpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.gcpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")



  vec.tx.ctpy <- unname(colMeans(52 * sim$epi$txCT / (sim$epi$num * sim$epi$prev.ct), na.rm = TRUE))
  df$cttxpy[i] <- paste0(round(quantile(vec.tx.ctpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.tx.ctpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.tx.ctpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
  vec.tx.ctpy.g1 <- unname(colMeans(52 * sim$epi$txCT.tttraj1 / (sim$epi$tt.traj.ct1 * sim$epi$prev.ct.tttraj1)))
  df.prev.ct.tttraj2 <- sim$epi$txCT.tttraj2[1:521,] / (sim$epi$tt.traj.ct2[1:521,] * sim$epi$prev.ct.tttraj2[1:521,])

  # Remove NaNs
  for (k in 1:ncol(df.prev.ct.tttraj2)) {
    df.prev.ct.tttraj2[which(is.nan(df.prev.ct.tttraj2[, k])), k] <- 0.0
  }
  vec.tx.ctpy.g2 <- unname(colMeans(52 * df.prev.ct.tttraj2))
  df$cttxpy.g1[i] <- paste0(round(quantile(vec.tx.ctpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.ctpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.ctpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$cttxpy.g2[i] <- paste0(round(quantile(vec.tx.ctpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.ctpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.ctpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  #STI
  vec.tx.stipy <- unname(colMeans(52 * (sim$epi$txGC + sim$epi$txCT) /
                                    (sim$epi$num * sim$epi$prev.sti)))
  vec.tx.stipy.g1 <- unname(colMeans(52 * (sim$epi$txGC.tttraj1 + sim$epi$txCT.tttraj1) /
                                       (sim$epi$tt.traj.sti1 * sim$epi$prev.sti.tttraj1)))
  df.prev.sti.tttraj2 <- (sim$epi$txGC.tttraj2[1:521,] + sim$epi$txCT.tttraj2[1:521, ]) / (sim$epi$tt.traj.sti2[1:521,] * sim$epi$prev.sti.tttraj2[1:521,])

  # Remove NaNs
  for (j in 1:ncol(df.prev.sti.tttraj2)) {

    df.prev.sti.tttraj2[which(is.nan(df.prev.sti.tttraj2[, j])), j] <- 0.0

  }
  vec.tx.stipy.g2 <- unname(colMeans(52 * df.prev.sti.tttraj2))

  df$txperpy.sti[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  df$txperpy.sti.g1[i] <- paste0(round(quantile(vec.tx.stipy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                 " (", round(quantile(vec.tx.stipy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                 " - ", round(quantile(vec.tx.stipy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                 ")")

  df$txperpy.sti.g2[i] <- paste0(round(quantile(vec.tx.stipy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                 " (", round(quantile(vec.tx.stipy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                 " - ", round(quantile(vec.tx.stipy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                 ")")
  # Number needed to treat
  # incremental tests / number of infections averted
  vec.gc.nnt <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - incid.gc)
  vec.gc.nnt.g1 <- (gc.asympt.tests.g1 - tests.gc.base.g1) / (incid.base.gc.g1 - incid.gc.g1)
  vec.gc.nnt.g2 <- (gc.asympt.tests.g2 - tests.gc.base.g2) / (incid.base.gc.g2 - incid.gc.g2)

  vec.ct.nnt <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - incid.ct)
  vec.ct.nnt.g1 <- (ct.asympt.tests.g1 - tests.ct.base.g1) / (incid.base.ct.g1 - incid.ct.g1)
  vec.ct.nnt.g2 <- (ct.asympt.tests.g2 - tests.ct.base.g2) / (incid.base.ct.g2 - incid.ct.g2)

  vec.gcct.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base.gcct - incid.gcct)
  vec.gcct.nnt.g1 <- (sti.asympt.tests.g1  - tests.sti.base.g1) / (incid.base.gcct.g1 - incid.gcct.g1)
  vec.gcct.nnt.g2 <- (sti.asympt.tests.g2  - tests.sti.base.g2) / (incid.base.gcct.g2 - incid.gcct.g2)

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

  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")
  df$gcct.nnt.g1[i] <- paste0(round(quantile(vec.gcct.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.gcct.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.gcct.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
  df$gcct.nnt.g2[i] <- paste0(round(quantile(vec.gcct.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.gcct.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.gcct.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
  cat("*")

}

df

write.csv(df, "analysis/STD Table S13.csv")
