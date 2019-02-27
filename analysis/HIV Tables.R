## HIV Table 1 ---------------------------------------------------------

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

hiv.paf.base <- as.numeric(quantile(unname(colMeans(tail(sim.base$epi$ir100, 52), na.rm = TRUE)), 0.5))

incid.base.gc <- unname(colSums(sim.base$epi$incid.rgc)) + unname(colSums(sim.base$epi$incid.ugc))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))

incid.base.ct <- unname(colSums(sim.base$epi$incid.rct)) + unname(colSums(sim.base$epi$incid.uct))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))

incid.base.gcct <- incid.base.ct + incid.base.gc
tests.sti.base <- tests.gc.base + tests.ct.base

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
gc.nnt <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.asympt.tests.py <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
sti.asympt.tests.py <- rep(NA, length(sims))
txperpy.sti <- rep(NA, length(sims))
gcct.nnt <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.pia <- rep(NA, length(sims))
hiv.tests <- rep(NA, length(sims))
hiv.nnt <- rep(NA, length(sims))
hiv.paf <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, partcut,

                 # HIV
                 hiv.incid, hiv.pia, hiv.nnt, hiv.tests, hiv.paf,
                 # STI
                 gc.incid, gc.pia, gc.nnt, gc.asympt.tests.py,
                 ct.incid, ct.pia, ct.nnt, ct.asympt.tests.py,
                 gcct.incid, gcct.pia, gcct.nnt, txperpy.sti, sti.asympt.tests.py

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
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52), na.rm = TRUE))
  df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52), na.rm = TRUE))
  df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  vec.ir.gcct <- vec.ir.gc + vec.ir.ct
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")

  # PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.rgc)) + unname(colSums(sim$epi$incid.ugc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.rct)) + unname(colSums(sim$epi$incid.uct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.gcct <- incid.gc + incid.ct
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")


  #STI Tx per year
  vec.tx.stipy <- unname(colMeans(52 * (sim$epi$txGC + sim$epi$txCT) /
                                    (sim$epi$num * sim$epi$prev.sti)))
  df$txperpy.sti[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  #Tests (Cumulative or over first year?)
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  gc.sympt.tests <- unname(colSums(sim$epi$GCsympttests, na.rm = TRUE))
  gc.tests <- gc.asympt.tests + gc.sympt.tests

  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  ct.sympt.tests <- unname(colSums(sim$epi$CTsympttests, na.rm = TRUE))
  ct.tests <- ct.asympt.tests + ct.sympt.tests

  sti.asympt.tests <- gc.asympt.tests + ct.asympt.tests
  sti.sympt.tests <- gc.sympt.tests + ct.sympt.tests
  sti.tests <- sti.asympt.tests + sti.sympt.tests

  py <- unname(colSums(sim$epi$num, na.rm = TRUE))

  gc.asympt.tests.py <-  52 * (gc.asympt.tests / py)
  ct.asympt.tests.py <-  52 * (ct.asympt.tests / py)
  sti.asympt.tests.py <-  52 * (sti.asympt.tests / py)

  df$gc.asympt.tests.py[i] <- paste0(round(quantile(gc.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(gc.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(gc.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")

  df$ct.asympt.tests.py[i] <- paste0(round(quantile(ct.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(ct.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(ct.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")
  df$sti.asympt.tests.py[i] <- paste0(round(quantile(sti.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                      " (", round(quantile(sti.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                      " - ", round(quantile(sti.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                      ")")

  # NNT
  vec.gc.nnt <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - incid.gc)
  vec.ct.nnt <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - incid.ct)
  vec.gcct.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base.gcct - incid.gcct)
  df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  # HIV incidence over the last year
  vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 52), na.rm = TRUE))
  df$hiv.incid[i] <- paste0(round(quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")

  # PAF
  df$hiv.paf[i] <- paste0(round(100 * ((hiv.paf.base -
                                   as.numeric(quantile(vec.ir.hiv, 0.5))) /
    (as.numeric(quantile(vec.ir.hiv, 0.5)))), 2), " %")

  # HIV PIA (Cumulative)
  incid.hiv <- unname(colSums(sim$epi$incid))
  vec.nia.hiv <- incid.base - incid.hiv
  vec.pia.hiv <- vec.nia.hiv/incid.base

  df$hiv.pia[i] <- paste0(round(quantile(vec.pia.hiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  # NNT
  vec.hiv.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base - incid.hiv)
  df$hiv.nnt[i] <- paste0(round(quantile(vec.hiv.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.hiv.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.hiv.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  # HIV Tests
  hivtests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  df$hiv.tests[i] <- paste0(round(quantile(hivtests, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(hivtests, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(hivtests, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  cat("*")

}

df

write.csv(df, "analysis/HIV Table 1.csv")

## HIV Table 2 ---------------------------------------------------------

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

hiv.paf.base <- as.numeric(quantile(unname(colMeans(tail(sim.base$epi$ir100, 52), na.rm = TRUE)), 0.5))

incid.base.gc <- unname(colSums(sim.base$epi$incid.rgc)) + unname(colSums(sim.base$epi$incid.ugc))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))

incid.base.ct <- unname(colSums(sim.base$epi$incid.rct)) + unname(colSums(sim.base$epi$incid.uct))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))

incid.base.gcct <- incid.base.ct + incid.base.gc
tests.sti.base <- tests.gc.base + tests.ct.base

# Varying Lower-Risk Coverage
# 9002, 9004, 9006, 9008
# : Annual = 10%, 20%, 30%, 40% increase, 364 days, HR = 0%, 182 days
# Varying Higher-Risk Coverage
# 9010, 9012, 9014, 9016, 9018, 9020, 9022, 9024, 9026, 9028
#  Higher-risk = 0.1 - 1.0 by 0.1, 182 days, 364 days

# Newer way:
sims <- c(9000, 9002, 9004, 9006, 9008,
          9010, 9012, 9014, 9016, 9018, 9020, 9022, 9024, 9026, 9028)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
gc.asympt.tests.py <- rep(NA, length(sims))
gc.nnt <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.asympt.tests.py <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
sti.asympt.tests.py <- rep(NA, length(sims))
txperpy.sti <- rep(NA, length(sims))
gcct.nnt <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.pia <- rep(NA, length(sims))
hiv.tests <- rep(NA, length(sims))
hiv.nnt <- rep(NA, length(sims))
hiv.paf <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint,

                 # HIV
                 hiv.incid, hiv.pia, hiv.nnt, hiv.tests, hiv.paf,

                 # STI
                 gc.incid, gc.pia, gc.nnt, gc.asympt.tests.py,
                 ct.incid, ct.pia, ct.nnt, ct.asympt.tests.py,
                 gcct.incid, gcct.pia, gcct.nnt, txperpy.sti, sti.asympt.tests.py

)

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int

  # Incidence Rate over last year
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52), na.rm = TRUE))
  df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52), na.rm = TRUE))
  df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  vec.ir.gcct <- vec.ir.gc + vec.ir.ct
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")

  # PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.rgc)) + unname(colSums(sim$epi$incid.ugc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.rct)) + unname(colSums(sim$epi$incid.uct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.gcct <- incid.gc + incid.ct
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")


  #Tests (Cumulative or over first year?)
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  gc.sympt.tests <- unname(colSums(sim$epi$GCsympttests, na.rm = TRUE))
  gc.tests <- gc.asympt.tests + gc.sympt.tests

  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  ct.sympt.tests <- unname(colSums(sim$epi$CTsympttests, na.rm = TRUE))
  ct.tests <- ct.asympt.tests + ct.sympt.tests

  sti.asympt.tests <- gc.asympt.tests + ct.asympt.tests
  sti.sympt.tests <- gc.sympt.tests + ct.sympt.tests
  sti.tests <- sti.asympt.tests + sti.sympt.tests

  py <- unname(colSums(sim$epi$num, na.rm = TRUE))

  gc.asympt.tests.py <-  52 * (gc.asympt.tests / py)
  ct.asympt.tests.py <-  52 * (ct.asympt.tests / py)
  sti.asympt.tests.py <-  52 * (sti.asympt.tests / py)

  df$gc.asympt.tests.py[i] <- paste0(round(quantile(gc.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(gc.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(gc.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")

  df$ct.asympt.tests.py[i] <- paste0(round(quantile(ct.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(ct.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(ct.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")
  df$sti.asympt.tests.py[i] <- paste0(round(quantile(sti.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                      " (", round(quantile(sti.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                      " - ", round(quantile(sti.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                      ")")

  # NNT
  vec.gc.nnt <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - incid.gc)
  vec.ct.nnt <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - incid.ct)
  vec.gcct.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base.gcct - incid.gcct)
  df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  #STI Tx per year
  vec.tx.stipy <- unname(colMeans(52 * (sim$epi$txGC + sim$epi$txCT) /
                                    (sim$epi$num * sim$epi$prev.sti)))
  df$txperpy.sti[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  # HIV incidence over the last year
  vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 52), na.rm = TRUE))
  df$hiv.incid[i] <- paste0(round(quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  # PAF
  df$hiv.paf[i] <- paste0(round(100 * ((hiv.paf.base -
                                          as.numeric(quantile(vec.ir.hiv, 0.5))) /
                                         (as.numeric(quantile(vec.ir.hiv, 0.5)))), 2), " %")

  # HIV PIA (Cumulative)
  incid.hiv <- unname(colSums(sim$epi$incid))
  vec.nia.hiv <- incid.base - incid.hiv
  vec.pia.hiv <- vec.nia.hiv/incid.base

  df$hiv.pia[i] <- paste0(round(quantile(vec.pia.hiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

  # NNT
  vec.hiv.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base - incid.hiv)
  df$hiv.nnt[i] <- paste0(round(quantile(vec.hiv.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.hiv.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.hiv.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

  # HIV Tests
  hivtests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  df$hiv.tests[i] <- paste0(round(quantile(hivtests, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(hivtests, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(hivtests, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  cat("*")

}

df

write.csv(df, "analysis/HIV Table 2.csv")

## HIV Table 3 ---------------------------------------------------------

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

hiv.paf.base <- as.numeric(quantile(unname(colMeans(tail(sim.base$epi$ir100, 52), na.rm = TRUE)), 0.5))

incid.base.gc <- unname(colSums(sim.base$epi$incid.rgc)) + unname(colSums(sim.base$epi$incid.ugc))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))

incid.base.ct <- unname(colSums(sim.base$epi$incid.rct)) + unname(colSums(sim.base$epi$incid.uct))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))

incid.base.gcct <- incid.base.ct + incid.base.gc
tests.sti.base <- tests.gc.base + tests.ct.base

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
gc.asympt.tests.py <- rep(NA, length(sims))
gc.nnt <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.asympt.tests.py <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
sti.asympt.tests.py <- rep(NA, length(sims))
txperpy.sti <- rep(NA, length(sims))
gcct.nnt <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.pia <- rep(NA, length(sims))
hiv.tests <- rep(NA, length(sims))
hiv.nnt <- rep(NA, length(sims))
hiv.paf <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, txprob,

                 # HIV
                 hiv.incid, hiv.pia, hiv.nnt, hiv.tests, hiv.paf,

                 # STI
                 gc.incid, gc.pia, gc.nnt, gc.asympt.tests.py,
                 ct.incid, ct.pia, ct.nnt, ct.asympt.tests.py,
                 gcct.incid, gcct.pia, gcct.nnt, txperpy.sti, sti.asympt.tests.py

)

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int
  df$txprob[i] <- sim$param$gc.asympt.prob.tx

  # Incidence Rate over last year
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52), na.rm = TRUE))
  df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52), na.rm = TRUE))
  df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  vec.ir.gcct <- vec.ir.gc + vec.ir.ct
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")

  # PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.rgc)) + unname(colSums(sim$epi$incid.ugc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.rct)) + unname(colSums(sim$epi$incid.uct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.gcct <- incid.gc + incid.ct
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")


  #Tests (Cumulative or over first year?)
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  gc.sympt.tests <- unname(colSums(sim$epi$GCsympttests, na.rm = TRUE))
  gc.tests <- gc.asympt.tests + gc.sympt.tests

  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  ct.sympt.tests <- unname(colSums(sim$epi$CTsympttests, na.rm = TRUE))
  ct.tests <- ct.asympt.tests + ct.sympt.tests

  sti.asympt.tests <- gc.asympt.tests + ct.asympt.tests
  sti.sympt.tests <- gc.sympt.tests + ct.sympt.tests
  sti.tests <- sti.asympt.tests + sti.sympt.tests

  py <- unname(colSums(sim$epi$num, na.rm = TRUE))

  gc.asympt.tests.py <-  52 * (gc.asympt.tests / py)
  ct.asympt.tests.py <-  52 * (ct.asympt.tests / py)
  sti.asympt.tests.py <-  52 * (sti.asympt.tests / py)

  df$gc.asympt.tests.py[i] <- paste0(round(quantile(gc.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(gc.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(gc.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")

  df$ct.asympt.tests.py[i] <- paste0(round(quantile(ct.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(ct.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(ct.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")
  df$sti.asympt.tests.py[i] <- paste0(round(quantile(sti.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                      " (", round(quantile(sti.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                      " - ", round(quantile(sti.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                      ")")

  # NNT
  vec.gc.nnt <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - incid.gc)
  vec.ct.nnt <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - incid.ct)
  vec.gcct.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base.gcct - incid.gcct)
  df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  #STI Tx per year
  vec.tx.stipy <- unname(colMeans(52 * (sim$epi$txGC + sim$epi$txCT) /
                                    (sim$epi$num * sim$epi$prev.sti)))
  df$txperpy.sti[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  # HIV incidence over the last year
  vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 52), na.rm = TRUE))
  df$hiv.incid[i] <- paste0(round(quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  # PAF
  df$hiv.paf[i] <- paste0(round(100 * ((hiv.paf.base -
                                          as.numeric(quantile(vec.ir.hiv, 0.5))) /
                                         (as.numeric(quantile(vec.ir.hiv, 0.5)))), 2), " %")

  # HIV PIA (Cumulative)
  incid.hiv <- unname(colSums(sim$epi$incid))
  vec.nia.hiv <- incid.base - incid.hiv
  vec.pia.hiv <- vec.nia.hiv/incid.base

  df$hiv.pia[i] <- paste0(round(quantile(vec.pia.hiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

  # NNT
  vec.hiv.nnt <- (sti.asympt.tests  - tests.sti.base) / (incid.base - incid.hiv)
  df$hiv.nnt[i] <- paste0(round(quantile(vec.hiv.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.hiv.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.hiv.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

  # HIV Tests
  hivtests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  df$hiv.tests[i] <- paste0(round(quantile(hivtests, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(hivtests, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(hivtests, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
  cat("*")

}

df

write.csv(df, "analysis/HIV Table 3.csv")
