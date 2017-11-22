## STI Testing Guidelines New Table 1 - NG/CT and Syphilis
# Sensitivity Analyses

## Compare to baseline ---------------------------------------------------------

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
# Reference scenario here
load("data/sim.n3000.rda")
#load("data/sim.n3191.rda")
sim.base <- sim

incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct))
incid.base.gcct.g1 <- unname(colSums(sim.base$epi$incid.gcct.tttraj1))
incid.base.gcct.g2 <- unname(colSums(sim.base$epi$incid.gcct.tttraj2))
tests.gcct.base <- unname(colSums(sim.base$epi$GCasympttests)) + unname(colSums(sim.base$epi$CTasympttests))
tests.gcct.base.g1 <- unname(colSums(sim.base$epi$GCasympttests.tttraj1)) + unname(colSums(sim.base$epi$CTasympttests.tttraj1))
tests.gcct.base.g2 <- unname(colSums(sim.base$epi$GCasympttests.tttraj2)) + unname(colSums(sim.base$epi$CTasympttests.tttraj2))

incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
incid.base.syph.g1 <- unname(colSums(sim.base$epi$incid.syph.tttraj1))
incid.base.syph.g2 <- unname(colSums(sim.base$epi$incid.syph.tttraj2))
tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests))
tests.syph.base.g1 <- unname(colSums(sim.base$epi$syphasympttests.tttraj1))
tests.syph.base.g2 <- unname(colSums(sim.base$epi$syphasympttests.tttraj2))

# Baseline compared to 5% HR coverage

# Screening Intervals:
# 3191 as reference
# 3189, 3190, 3191, 3192, 3193: Annual = 182 days, 273 days, 364 days (ref), 448 days, 539 days, HR = 5%, Ann = baseline,
# 3194-3198: Higher-risk = 28 days, 91 days, 182 days (ref), 273 days, 364 days, HR = 5%, Ann = Baseline, 364 days
# 3191 and 3196 should be approx equal
# Treatment Completion - for appendix
# 3199, 3204, 3209, 3214: Annual = 0.0 - 1.0 (ref) by 0.25, 364 days, HR = 0%, 182 days
#
# Partner Cutoff for Higher-Risk
# 3221:3229 Higher-risk = 1 (ref) to 10 by 1

# Newer way:
sims <- c(3000, 3189:3193, 3194:3195, 3191, 3197:3198, 3191, 3221:3229)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
partcut <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
gcct.asympt.tests.py <- rep(NA, length(sims))
gcct.asympt.tests <- rep(NA, length(sims))
# gccttx <- rep(NA, length(sims))
# gccttxpy <- rep(NA, length(sims))

gcct.incid.g1 <- rep(NA, length(sims))
gcct.pia.g1 <- rep(NA, length(sims))
gcct.asympt.tests.py.g1 <- rep(NA, length(sims))
gcct.asympt.tests.g1 <- rep(NA, length(sims))
# gccttx.g1 <- rep(NA, length(sims))
# gccttxpy.g1 <- rep(NA, length(sims))

gcct.incid.g2 <- rep(NA, length(sims))
gcct.pia.g2 <- rep(NA, length(sims))
gcct.asympt.tests.py.g2 <- rep(NA, length(sims))
gcct.asympt.tests.g2 <- rep(NA, length(sims))
# gccttx.g2 <- rep(NA, length(sims))
# gccttxpy.g2 <- rep(NA, length(sims))

gcct.nnt <- rep(NA, length(sims))
gcct.nnt.g1 <- rep(NA, length(sims))
gcct.nnt.g2 <- rep(NA, length(sims))

syph.incid <- rep(NA, length(sims))
syph.pia <- rep(NA, length(sims))
syph.asympt.tests.py <- rep(NA, length(sims))
syph.asympt.tests <- rep(NA, length(sims))
# syphtx <- rep(NA, length(sims))
# syphearlytx <- rep(NA, length(sims))
# syphlatetx <- rep(NA, length(sims))
# syphtxpy <- rep(NA, length(sims))
# syphearlytxpy <- rep(NA, length(sims))
# syphlatetxpy <- rep(NA, length(sims))

syph.incid.g1 <- rep(NA, length(sims))
syph.pia.g1 <- rep(NA, length(sims))
syph.asympt.tests.g1 <- rep(NA, length(sims))
syph.asympt.tests.py.g1 <- rep(NA, length(sims))
# syphtx.g1 <- rep(NA, length(sims))
# syphearlytx.g1 <- rep(NA, length(sims))
# syphlatetx.g1 <- rep(NA, length(sims))
# syphtxpy.g1 <- rep(NA, length(sims))
# syphearlytxpy.g1 <- rep(NA, length(sims))
# syphlatetxpy.g1 <- rep(NA, length(sims))

syph.incid.g2 <- rep(NA, length(sims))
syph.pia.g2 <- rep(NA, length(sims))
syph.asympt.tests.g2 <- rep(NA, length(sims))
syph.asympt.tests.py.g2 <- rep(NA, length(sims))
# syphtx.g2 <- rep(NA, length(sims))
# syphearlytx.g2 <- rep(NA, length(sims))
# syphlatetx.g2 <- rep(NA, length(sims))
# syphtxpy.g2 <- rep(NA, length(sims))
# syphearlytxpy.g2 <- rep(NA, length(sims))
# syphlatetxpy.g2 <- rep(NA, length(sims))

syph.nnt <- rep(NA, length(sims))
syph.nnt.g1 <- rep(NA, length(sims))
syph.nnt.g2 <- rep(NA, length(sims))


df <- data.frame(anncov, hrcov, annint, hrint, partcut,

                 # Overall
                 gcct.incid, gcct.pia, gcct.nnt, gcct.asympt.tests.py, gcct.asympt.tests,# gccttx, gccttxpy,
                 syph.incid, syph.pia, syph.nnt, syph.asympt.tests.py, syph.asympt.tests,

                 # Group 1
                 gcct.incid.g1, gcct.pia.g1, gcct.nnt.g1,
                 gcct.asympt.tests.py.g1, gcct.asympt.tests.g1, #gccttx.g1, gccttxpy.g1,
                 syph.incid.g1, syph.pia.g1, syph.nnt.g1,
                 syph.asympt.tests.py.g1, syph.asympt.tests.g1, #syphtx.g1, syphearlytx.g1,
                 #syphlatetx.g1, syphtxpy.g1, syphearlytxpy.g1, syphlatetxpy.g1,

                 # Group 2
                 gcct.incid.g2, gcct.pia.g2, gcct.nnt.g2,
                 gcct.asympt.tests.py.g2, gcct.asympt.tests.g2, gccttx.g2, gccttxpy.g2,
                 syph.incid.g2, syph.pia.g2, syph.nnt.g2,
                 syph.asympt.tests.py.g2, syph.asympt.tests.g2, #syphtx.g2,
                # syphearlytx.g2, syphlatetx.g2, syphtxpy.g2, syphearlytxpy.g2, syphlatetxpy.g2
)

for (i in seq_along(sims)) {

  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int
  df$partcut[i] <- sim$param$partnercutoff

  # Incidence Rate over last year
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gcct, 52)))
  vec.ir.gc.g1 <- unname(colMeans(tail(sim$epi$ir100.gcct.tttraj1, 52)))
  vec.ir.gc.g2 <- unname(colMeans(tail(sim$epi$ir100.gcct.tttraj2, 52)))
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$gcct.incid.g1[i] <- paste0(round(quantile(vec.ir.gc.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.gc.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.gc.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$gcct.incid.g2[i] <- paste0(round(quantile(vec.ir.gc.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.gc.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.gc.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  vec.ir.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
  vec.ir.syph.g1 <- unname(colMeans(tail(sim$epi$ir100.syph.tttraj1, 52)))
  vec.ir.syph.g2 <- unname(colMeans(tail(sim$epi$ir100.syph.tttraj2, 52)))
  df$syph.incid[i] <- paste0(round(quantile(vec.ir.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")
  df$syph.incid.g1[i] <- paste0(round(quantile(vec.ir.syph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                " (", round(quantile(vec.ir.syph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                " - ", round(quantile(vec.ir.syph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                ")")
  df$syph.incid.g2[i] <- paste0(round(quantile(vec.ir.syph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                " (", round(quantile(vec.ir.syph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                " - ", round(quantile(vec.ir.syph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                ")")

  # PIA (Cumulative)
  incid.gcct <- unname(colSums(sim$epi$incid.gcct))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  incid.gcct.g1 <- unname(colSums(sim$epi$incid.gcct.tttraj1))
  vec.nia.gcct.g1 <- round(incid.base.gcct.g1 - incid.gcct.g1, 1)
  vec.pia.gcct.g1 <- vec.nia.gcct.g1/incid.gcct.g1

  incid.gcct.g2 <- unname(colSums(sim$epi$incid.gcct.tttraj2))
  vec.nia.gcct.g2 <- round(incid.base.gcct.g2 - incid.gcct.g2, 1)
  vec.pia.gcct.g2 <- vec.nia.gcct.g2/incid.gcct.g2

  incid.syph <- unname(colSums(sim$epi$incid.syph))
  vec.nia.syph <- incid.base.syph - incid.syph
  vec.pia.syph <- vec.nia.syph/incid.base.syph

  incid.syph.g1 <- unname(colSums(sim$epi$incid.syph.tttraj1))
  vec.nia.syph.g1 <- round(incid.base.syph.g1 - incid.syph.g1, 1)
  vec.pia.syph.g1 <- vec.nia.syph.g1/incid.syph.g1

  incid.syph.g2 <- unname(colSums(sim$epi$incid.syph.tttraj2))
  vec.nia.syph.g2 <- round(incid.base.syph.g2 - incid.syph.g2, 1)
  vec.pia.syph.g2 <- vec.nia.syph.g1/incid.syph.g2

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

  df$syph.pia[i] <- paste0(round(quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  df$syph.pia.g1[i] <- paste0(round(quantile(vec.pia.syph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.pia.syph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.pia.syph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$syph.pia.g2[i] <- paste0(round(quantile(vec.pia.syph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.pia.syph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.pia.syph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  #Tests (Cumulative or over first year?)
  gcct.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE)) + unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  gcct.sympt.tests <- unname(colSums(sim$epi$GCsympttests, na.rm = TRUE)) + unname(colSums(sim$epi$CTsympttests, na.rm = TRUE))
  gcct.tests <- gcct.asympt.tests + gcct.sympt.tests

  gcct.asympt.tests.g1 <- unname(colSums(sim$epi$GCasympttests.tttraj1, na.rm = TRUE)) + unname(colSums(sim$epi$CTasympttests.tttraj1, na.rm = TRUE))
  gcct.sympt.tests.g1 <- unname(colSums(sim$epi$GCsympttests.tttraj1, na.rm = TRUE)) + unname(colSums(sim$epi$CTsympttests.tttraj1, na.rm = TRUE))
  gcct.tests.g1 <- gcct.asympt.tests.g1 + gcct.sympt.tests.g1

  gcct.asympt.tests.g2 <- unname(colSums(sim$epi$GCasympttests.tttraj2, na.rm = TRUE)) + unname(colSums(sim$epi$CTasympttests.tttraj2, na.rm = TRUE))
  gcct.sympt.tests.g2 <- unname(colSums(sim$epi$GCsympttests.tttraj2, na.rm = TRUE)) + unname(colSums(sim$epi$CTsympttests.tttraj2, na.rm = TRUE))
  gcct.tests.g2 <- gcct.asympt.tests.g2 + gcct.sympt.tests.g2

  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
  syph.sympt.tests <- unname(colSums(sim$epi$syphsympttests, na.rm = TRUE))
  syph.tests <- syph.asympt.tests + syph.sympt.tests

  syph.asympt.tests.g1 <- unname(colSums(sim$epi$syphasympttests.tttraj1, na.rm = TRUE))
  syph.sympt.tests.g1 <- unname(colSums(sim$epi$syphsympttests.tttraj1, na.rm = TRUE))
  syph.tests.g1 <- syph.asympt.tests.g1 + syph.sympt.tests.g1

  syph.asympt.tests.g2 <- unname(colSums(sim$epi$syphasympttests.tttraj2, na.rm = TRUE))
  syph.sympt.tests.g2 <- unname(colSums(sim$epi$syphsympttests.tttraj2, na.rm = TRUE))
  syph.tests.g2 <- syph.asympt.tests.g2 + syph.sympt.tests.g2

  py <- unname(colSums(sim$epi$num, na.rm = TRUE))
  py.g1 <- unname(colSums(sim$epi$tt.traj.sti1, na.rm = TRUE))
  py.g2 <- unname(colSums(sim$epi$tt.traj.sti2, na.rm = TRUE))

  gc.asympt.tests.py <-  52 * (gcct.asympt.tests / py)
  gc.asympt.tests.py.g1 <-  52 * (gcct.asympt.tests.g1 / py.g1)
  gc.asympt.tests.py.g2 <-  52 * (gcct.asympt.tests.g2 / py.g2)
  syph.asympt.tests.py <-  52 * (syph.asympt.tests / py)
  syph.asympt.tests.py.g1 <-  52 * (syph.asympt.tests.g1 / py.g1)
  syph.asympt.tests.py.g2 <-  52 * (syph.asympt.tests.g2 / py.g2)

  df$gcct.asympt.tests[i] <- paste0(round(quantile(gcct.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                  " (", round(quantile(gcct.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                  " - ", round(quantile(gcct.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                  ")")
  df$gcct.asympt.tests.g1[i] <- paste0(round(quantile(gcct.asympt.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                     " (", round(quantile(gcct.asympt.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                     " - ", round(quantile(gcct.asympt.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                     ")")
  df$gcct.asympt.tests.g2[i] <- paste0(round(quantile(gcct.asympt.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                     " (", round(quantile(gcct.asympt.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                     " - ", round(quantile(gcct.asympt.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                     ")")
  df$gcct.asympt.tests.py[i] <- paste0(round(quantile(gcct.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                     " (", round(quantile(gcct.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                     " - ", round(quantile(gcct.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                     ")")

  df$gcct.asympt.tests.py.g1[i] <- paste0(round(quantile(gcct.asympt.tests.py.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(gcct.asympt.tests.py.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(gcct.asympt.tests.py.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")
  df$gcct.asympt.tests.py.g2[i] <- paste0(round(quantile(gcct.asympt.tests.py.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(gcct.asympt.tests.py.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(gcct.asympt.tests.py.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")

  df$syph.asympt.tests[i] <- paste0(round(quantile(syph.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                    " (", round(quantile(syph.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                    " - ", round(quantile(syph.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                    ")")
  df$syph.asympt.tests.g1[i] <- paste0(round(quantile(syph.asympt.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                       " (", round(quantile(syph.asympt.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                       " - ", round(quantile(syph.asympt.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                       ")")
  df$syph.asympt.tests.g2[i] <- paste0(round(quantile(syph.asympt.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                                       " (", round(quantile(syph.asympt.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                                       " - ", round(quantile(syph.asympt.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                                       ")")
  df$syph.asympt.tests.py[i] <- paste0(round(quantile(syph.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                       " (", round(quantile(syph.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                       " - ", round(quantile(syph.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                       ")")

  df$syph.asympt.tests.py.g1[i] <- paste0(round(quantile(syph.asympt.tests.py.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                          " (", round(quantile(syph.asympt.tests.py.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                          " - ", round(quantile(syph.asympt.tests.py.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                          ")")
  df$syph.asympt.tests.py.g2[i] <- paste0(round(quantile(syph.asympt.tests.py.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                          " (", round(quantile(syph.asympt.tests.py.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                          " - ", round(quantile(syph.asympt.tests.py.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                          ")")

  # Number needed to treat
  vec.gc.nnt <- (gcct.asympt.tests - tests.gcct.base) / (incid.base.gcct - incid.gcct)
  vec.gc.nnt.g1 <- (gcct.asympt.tests.g1 - tests.gcct.base.g1) / (incid.base.gcct.g1 - incid.gcct.g1)
  vec.gc.nnt.g2 <- (gcct.asympt.tests.g2 - tests.gcct.base.g2) / (incid.base.gcct.g2 - incid.gcct.g2)

  vec.syph.nnt <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - incid.syph)
  vec.syph.nnt.g1 <- (syph.asympt.tests.g1  - tests.syph.base.g1) / (incid.base.syph.g1 - incid.syph.g1)
  vec.syph.nnt.g2 <- (syph.asympt.tests.g2  - tests.syph.base.g2) / (incid.base.syph.g2 - incid.syph.g2)

  df$gcct.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.nnt.g1[i] <- paste0(round(quantile(vec.gc.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.gc.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.gc.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")
  df$gcct.nnt.g2[i] <- paste0(round(quantile(vec.gc.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.gc.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.gc.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
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

  cat("*")

}

df

write.csv(df, "analysis/New STD Table 1.csv")

