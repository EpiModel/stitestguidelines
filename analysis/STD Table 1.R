## STI Testing Guidelines Table 2
# Sensitivity Analyses

## Compare to baseline ---------------------------------------------------------

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

# Base - No annual or high-risk
# Reference scenario here
load("data/followup/Guidelines Paper/sim.n9000.rda")
sim.base <- sim

incid.base <- unname(colSums(sim.base$epi$incid))
tests.base <- unname(colSums(sim.base$epi$hivtests.nprep))

incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
incid.base.gc.g1 <- unname(colSums(sim.base$epi$incid.gc.tttraj1))
incid.base.gc.g2 <- unname(colSums(sim.base$epi$incid.gc.tttraj2))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))
tests.gc.base.g1 <- unname(colSums(sim.base$epi$GCasympttests.tttraj1))
tests.gc.base.g2 <- unname(colSums(sim.base$epi$GCasympttests.tttraj2))

incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
incid.base.ct.g1 <- unname(colSums(sim.base$epi$incid.ct.tttraj1))
incid.base.ct.g2 <- unname(colSums(sim.base$epi$incid.ct.tttraj2))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))
tests.ct.base.g1 <- unname(colSums(sim.base$epi$CTasympttests.tttraj1))
tests.ct.base.g2 <- unname(colSums(sim.base$epi$CTasympttests.tttraj2))

incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct))
incid.base.gcct.g1 <- unname(colSums(sim.base$epi$incid.gcct.tttraj1))
incid.base.gcct.g2 <- unname(colSums(sim.base$epi$incid.gcct.tttraj2))
tests.sti.base <- unname(colSums(sim.base$epi$stiasympttests))
tests.sti.base.g1 <- unname(colSums(sim.base$epi$stiasympttests.tttraj1))
tests.sti.base.g2 <- unname(colSums(sim.base$epi$stiasympttests.tttraj2))

# Baseline compared to 5% HR coverage

# 9000 as baseline
# 9009 as reference (5% HR, 182 days, baseline SA, 364 days)
# 9029:9030, 9009, 9031:9032: Annual = 182 days, 273 days, 364 days, 448 days, 539 days, HR = 5%, Ann = baseline,
# 9033:9034, 9009, 9035:9036: Higher-risk = 28 days, 91 days, 182 days, 273 days, 364 days, HR = 5%, Ann = Baseline, 364 days
#
# Partner Cutoff for Higher-Risk
# 9009 - 1 partner, 5% hR
# 9040:9048 Higher-risk = 2 to 10 (ref = 9032)

# Newer way:
sims <- c(9000, 9029:9030, 9009, 9031:9032, 9033:9034, 9009, 9035:9036, 9009, 9037:9045)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
partcut <- rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
gc.asympt.tests.py <- rep(NA, length(sims))
gc.asympt.tests <- rep(NA, length(sims))
gctx <- rep(NA, length(sims))
gctxpy <- rep(NA, length(sims))

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

gc.nnt <- rep(NA, length(sims))
gc.nnt.g1 <- rep(NA, length(sims))
gc.nnt.g2 <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.asympt.tests.py <- rep(NA, length(sims))
ct.asympt.tests <- rep(NA, length(sims))
cttx <- rep(NA, length(sims))
cttxpy <- rep(NA, length(sims))

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

ct.nnt <- rep(NA, length(sims))
ct.nnt.g1 <- rep(NA, length(sims))
ct.nnt.g2 <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
gcct.pia.g1 <- rep(NA, length(sims))
gcct.pia.g2 <- rep(NA, length(sims))
sti.asympt.tests <- rep(NA, length(sims))
sti.asympt.tests.py <- rep(NA, length(sims))
tx.sti <- rep(NA, length(sims))
txperpy.sti <- rep(NA, length(sims))

gcct.incid.g1 <- rep(NA, length(sims))
sti.pia.g1 <- rep(NA, length(sims))
sti.asympt.tests.g1 <- rep(NA, length(sims))
sti.asympt.tests.py.g1 <- rep(NA, length(sims))
tx.sti.g1 <- rep(NA, length(sims))
txperpy.sti.g1 <- rep(NA, length(sims))

gcct.incid.g2 <- rep(NA, length(sims))
sti.pia.g2 <- rep(NA, length(sims))
sti.asympt.tests.g2 <- rep(NA, length(sims))
sti.asympt.tests.py.g2 <- rep(NA, length(sims))
tx.sti.g2 <- rep(NA, length(sims))
txperpy.sti.g2 <- rep(NA, length(sims))

gcct.nnt <- rep(NA, length(sims))
gcct.nnt.g1 <- rep(NA, length(sims))
gcct.nnt.g2 <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, partcut,

                 # Overall
                 gc.incid, gc.pia, gc.nnt, gc.asympt.tests.py, gc.asympt.tests, gctx, gctxpy,
                 ct.incid, ct.pia, ct.nnt, ct.asympt.tests.py, ct.asympt.tests, cttx, cttxpy,
                 gcct.incid, gcct.pia, gcct.nnt, txperpy.sti, sti.asympt.tests.py,
                 sti.asympt.tests, tx.sti,

                 # Group 1
                 gc.incid.g1, gc.pia.g1, gc.nnt.g1,
                 gc.asympt.tests.py.g1, gc.asympt.tests.g1, gctx.g1, gctxpy.g1,
                 ct.incid.g1, ct.pia.g1, ct.nnt.g1,
                 ct.asympt.tests.py.g1, ct.asympt.tests.g1, cttx.g1, cttxpy.g1,
                 gcct.incid.g2, gcct.pia.g1, gcct.nnt.g1, txperpy.sti.g1,
                 sti.asympt.tests.py.g1, sti.asympt.tests.g1, tx.sti.g1,

                 # Group 2
                 gc.incid.g2, gc.pia.g2, gc.nnt.g2,
                 gc.asympt.tests.py.g2, gc.asympt.tests.g2, gctx.g2, gctxpy.g2,
                 ct.incid.g2, ct.pia.g2, ct.nnt.g2,
                 ct.asympt.tests.py.g2, ct.asympt.tests.g2, cttx.g2, cttxpy.g2,

                 gcct.incid.g2, gcct.pia.g2, gcct.nnt.g2,  txperpy.sti.g2,
                 sti.asympt.tests.py.g2, sti.asympt.tests.g2, tx.sti.g2
                 )

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int
  df$partcut[i] <- sim$param$partnercutoff

  # Incidence Rate over last year
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

  vec.ir.gcct <- unname(colMeans(tail(sim$epi$ir100.gcct, 52)))
  vec.ir.gcct.g1 <- unname(colMeans(tail(sim$epi$ir100.gcct.tttraj1, 52)))
  vec.ir.gcct.g2 <- unname(colMeans(tail(sim$epi$ir100.gcct.tttraj2, 52)))
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

  # PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.gc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.gc.g1 <- unname(colSums(sim$epi$incid.gc.tttraj1))
  vec.nia.gc.g1 <- round(incid.base.gc.g1 - incid.gc.g1, 1)
  vec.pia.gc.g1 <- vec.nia.gc.g1/incid.gc.g1

  incid.gc.g2 <- unname(colSums(sim$epi$incid.gc.tttraj2))
  vec.nia.gc.g2 <- round(incid.base.gc.g2 - incid.gc.g2, 1)
  vec.pia.gc.g2 <- vec.nia.gc.g2/incid.gc.g2

  incid.ct <- unname(colSums(sim$epi$incid.ct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.ct.g1 <- unname(colSums(sim$epi$incid.ct.tttraj1))
  vec.nia.ct.g1 <- round(incid.base.ct.g1 - incid.ct.g1, 1)
  vec.pia.ct.g1 <- vec.nia.ct.g1/incid.ct.g1

  incid.ct.g2 <- unname(colSums(sim$epi$incid.ct.tttraj2))
  vec.nia.ct.g2 <- round(incid.base.ct.g2 - incid.ct.g2, 1)
  vec.pia.ct.g2 <- vec.nia.ct.g2/incid.ct.g2

  incid.gcct <- unname(colSums(sim$epi$incid.gcct))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  incid.gcct.g1 <- unname(colSums(sim$epi$incid.gcct.tttraj1))
  vec.nia.gcct.g1 <- round(incid.base.gcct.g1 - incid.gcct.g1, 1)
  vec.pia.gcct.g1 <- vec.nia.gcct.g1/incid.gcct.g1

  incid.gcct.g2 <- unname(colSums(sim$epi$incid.gcct.tttraj2))
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

  #Tests (Cumulative or over first year?)
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

  sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))
  sti.sympt.tests <- unname(colSums(sim$epi$stisympttests, na.rm = TRUE))
  sti.tests <- sti.asympt.tests + sti.sympt.tests

  sti.asympt.tests.g1 <- unname(colSums(sim$epi$stiasympttests.tttraj1, na.rm = TRUE))
  sti.sympt.tests.g1 <- unname(colSums(sim$epi$stisympttests.tttraj1, na.rm = TRUE))
  sti.tests.g1 <- sti.asympt.tests.g1 + sti.sympt.tests.g1

  sti.asympt.tests.g2 <- unname(colSums(sim$epi$stiasympttests.tttraj2, na.rm = TRUE))
  sti.sympt.tests.g2 <- unname(colSums(sim$epi$stisympttests.tttraj2, na.rm = TRUE))
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

  vec.tx.gc <- as.numeric(colSums(head(sim$epi$txGC, 52), na.rm = TRUE))
  vec.tx.gc.g1 <- as.numeric(colSums(head(sim$epi$txGC.tttraj1, 52), na.rm = TRUE))
  vec.tx.gc.g2 <- as.numeric(colSums(head(sim$epi$txGC.tttraj2, 52), na.rm = TRUE))

  vec.tx.ct <- as.numeric(colSums(head(sim$epi$txCT, 52), na.rm = TRUE))
  vec.tx.ct.g1 <- as.numeric(colSums(head(sim$epi$txCT.tttraj1, 52), na.rm = TRUE))
  vec.tx.ct.g2 <- as.numeric(colSums(head(sim$epi$txCT.tttraj2, 52), na.rm = TRUE))

  vec.tx.sti <- unname(colSums(sim$epi$txSTI))
  vec.tx.sti.g1 <- unname(colSums(sim$epi$txSTI.tttraj1))
  vec.tx.sti.g2 <- unname(colSums(sim$epi$txSTI.tttraj2))

  df$gctx[i] <- paste0(round(quantile(vec.tx.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                        " (", round(quantile(vec.tx.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                        " - ", round(quantile(vec.tx.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                        ")")
  df$gctx.g1[i] <- paste0(round(quantile(vec.tx.gc.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                       " (", round(quantile(vec.tx.gc.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                       " - ", round(quantile(vec.tx.gc.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                       ")")
  df$gctx.g2[i] <- paste0(round(quantile(vec.tx.gc.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                       " (", round(quantile(vec.tx.gc.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                       " - ", round(quantile(vec.tx.gc.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                       ")")

  df$cttx[i] <- paste0(round(quantile(vec.tx.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                       " (", round(quantile(vec.tx.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                       " - ", round(quantile(vec.tx.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                       ")")
  df$cttx.g1[i] <- paste0(round(quantile(vec.tx.ct.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                          " (", round(quantile(vec.tx.ct.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                          " - ", round(quantile(vec.tx.ct.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                          ")")
  df$cttx.g2[i] <- paste0(round(quantile(vec.tx.ct.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                          " (", round(quantile(vec.tx.ct.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                          " - ", round(quantile(vec.tx.ct.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                          ")")

  df$tx.sti[i] <- paste0(round(quantile(vec.tx.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(vec.tx.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(vec.tx.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
  df$tx.sti.g1[i] <- paste0(round(quantile(vec.tx.sti.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                     " (", round(quantile(vec.tx.sti.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                     " - ", round(quantile(vec.tx.sti.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                     ")")
  df$tx.sti.g2[i] <- paste0(round(quantile(vec.tx.sti.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                     " (", round(quantile(vec.tx.sti.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                     " - ", round(quantile(vec.tx.sti.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                     ")")

  #txSTI or txSTI_asympt? (Over first year or cumulative?)
  #Update n values to num.tttraj1 and num.tttraj2

  # GC
  vec.tx.gcpy <- unname(colMeans(52 * sim$epi$txGC / (sim$epi$num * sim$epi$prev.gc)))
  vec.tx.gcpy.g1 <- unname(colMeans(52 * sim$epi$txGC.tttraj1 / (sim$epi$tt.traj.gc1 * sim$epi$prev.gc.tttraj1)))
  df.prev.gc.tttraj2 <- sim$epi$txGC.tttraj2[1:521,] / (sim$epi$tt.traj.gc2[1:521,] * sim$epi$prev.gc.tttraj2[1:521,])

  # Remove NaNs
  for (j in 1:ncol(df.prev.gc.tttraj2)) {
    df.prev.gc.tttraj2[which(is.nan(df.prev.gc.tttraj2[, j])), j] <- 0.0
  }
  vec.tx.gcpy.g2 <- unname(colMeans(52 * df.prev.gc.tttraj2))

  # CT
  vec.tx.ctpy <- unname(colMeans(52 * sim$epi$txCT / (sim$epi$num * sim$epi$prev.ct)))
  vec.tx.ctpy.g1 <- unname(colMeans(52 * sim$epi$txCT.tttraj1 / (sim$epi$tt.traj.ct1 * sim$epi$prev.ct.tttraj1)))
  df.prev.ct.tttraj2 <- sim$epi$txCT.tttraj2[1:521,] / (sim$epi$tt.traj.ct2[1:521,] * sim$epi$prev.ct.tttraj2[1:521,])
  # Remove NaNs
  for (k in 1:ncol(df.prev.ct.tttraj2)) {
    df.prev.ct.tttraj2[which(is.nan(df.prev.ct.tttraj2[, k])), k] <- 0.0
  }
  vec.tx.ctpy.g2 <- unname(colMeans(52 * df.prev.ct.tttraj2))

  df$gctxpy[i] <- paste0(round(quantile(vec.tx.gcpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.tx.gcpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.tx.gcpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$gctxpy.g1[i] <- paste0(round(quantile(vec.tx.gcpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.gcpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.gcpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$gctxpy.g2[i] <- paste0(round(quantile(vec.tx.gcpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.gcpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.gcpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$cttxpy[i] <- paste0(round(quantile(vec.tx.ctpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.tx.ctpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.tx.ctpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$cttxpy.g1[i] <- paste0(round(quantile(vec.tx.ctpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.ctpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.ctpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$cttxpy.g2[i] <- paste0(round(quantile(vec.tx.ctpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.ctpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.ctpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  #STI
  vec.tx.stipy <- unname(colMeans(52 * sim$epi$txSTI / (sim$epi$num * sim$epi$prev.sti)))
  vec.tx.stipy.g1 <- unname(colMeans(52 * sim$epi$txSTI.tttraj1 / (sim$epi$tt.traj.sti1 * sim$epi$prev.sti.tttraj1)))
  df.prev.sti.tttraj2 <- sim$epi$txSTI.tttraj2[1:521,] / (sim$epi$tt.traj.sti2[1:521,] * sim$epi$prev.sti.tttraj2[1:521,])

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

write.csv(df, "analysis/STD Table 1.csv")

