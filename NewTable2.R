## STI Testing Guidelines New Table 2 - Syph and NG/CT
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

#incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct))
incid.base.gcct <- unname(colSums(sim.base$epi$incid.gc)) + unname(colSums(sim.base$epi$incid.ct))
#incid.base.gcct.g1 <- unname(colSums(sim.base$epi$incid.gcct.tttraj1))
incid.base.gcct.g1 <- unname(colSums(sim.base$epi$incid.gc.tttraj1)) + unname(colSums(sim.base$epi$incid.ct.tttraj1))
#incid.base.gcct.g2 <- unname(colSums(sim.base$epi$incid.gcct.tttraj2))
incid.base.gcct.g2 <- unname(colSums(sim.base$epi$incid.gcct.tttraj2)) + unname(colSums(sim.base$epi$incid.ct.tttraj2))
tests.gcct.base <- unname(colSums(sim.base$epi$GCasympttests)) + unname(colSums(sim.base$epi$CTasympttests))
tests.gcct.base.g1 <- unname(colSums(sim.base$epi$GCasympttests.tttraj1)) + unname(colSums(sim.base$epi$CTasympttests.tttraj1))
tests.gcct.base.g2 <- unname(colSums(sim.base$epi$GCasympttests.tttraj2)) + unname(colSums(sim.base$epi$CTasympttests.tttraj2))

incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
incid.base.syph.g1 <- unname(colSums(sim.base$epi$incid.syph.tttraj1))
incid.base.syph.g2 <- unname(colSums(sim.base$epi$incid.syph.tttraj2))
tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests))
tests.syph.base.g1 <- unname(colSums(sim.base$epi$syphasympttests.tttraj1))
tests.syph.base.g2 <- unname(colSums(sim.base$epi$syphasympttests.tttraj2))

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

tx.gcct.prop <- rep(NA, length(sims))
tx.gcct.propv2 <- rep(NA, length(sims))
gcct.infect.dur <- rep(NA, length(sims))
gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
gcct.asympt.tests.py <- rep(NA, length(sims))
gcct.asympt.tests <- rep(NA, length(sims))
gccttxpy <- rep(NA, length(sims))

gcct.incid.g1 <- rep(NA, length(sims))
gcct.pia.g1 <- rep(NA, length(sims))
gcct.asympt.tests.py.g1 <- rep(NA, length(sims))
gcct.asympt.tests.g1 <- rep(NA, length(sims))
gccttxpy.g1 <- rep(NA, length(sims))

gcct.incid.g2 <- rep(NA, length(sims))
gcct.pia.g2 <- rep(NA, length(sims))
gcct.asympt.tests.py.g2 <- rep(NA, length(sims))
gcct.asympt.tests.g2 <- rep(NA, length(sims))
gccttxpy.g2 <- rep(NA, length(sims))

gcct.nnt <- rep(NA, length(sims))
gcct.nnt.g1 <- rep(NA, length(sims))
gcct.nnt.g2 <- rep(NA, length(sims))

tx.syph.prop <- rep(NA, length(sims))
tx.syph.propv2 <- rep(NA, length(sims))
syph.infect.dur <- rep(NA, length(sims))
syph.incid <- rep(NA, length(sims))
syph.pia <- rep(NA, length(sims))
syph.asympt.tests.py <- rep(NA, length(sims))
syph.asympt.tests <- rep(NA, length(sims))
syphtxpy <- rep(NA, length(sims))
syphearlytxpy <- rep(NA, length(sims))
syphlatetxpy <- rep(NA, length(sims))

syph.incid.g1 <- rep(NA, length(sims))
syph.pia.g1 <- rep(NA, length(sims))
syph.asympt.tests.g1 <- rep(NA, length(sims))
syph.asympt.tests.py.g1 <- rep(NA, length(sims))
syphtxpy.g1 <- rep(NA, length(sims))
syphearlytxpy.g1 <- rep(NA, length(sims))
syphlatetxpy.g1 <- rep(NA, length(sims))

syph.incid.g2 <- rep(NA, length(sims))
syph.pia.g2 <- rep(NA, length(sims))
syph.asympt.tests.g2 <- rep(NA, length(sims))
syph.asympt.tests.py.g2 <- rep(NA, length(sims))
syphtxpy.g2 <- rep(NA, length(sims))
syphearlytxpy.g2 <- rep(NA, length(sims))
syphlatetxpy.g2 <- rep(NA, length(sims))

syph.nnt <- rep(NA, length(sims))
syph.nnt.g1 <- rep(NA, length(sims))
syph.nnt.g2 <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint,
                 # Overall
                 gcct.incid, gcct.pia, gcct.nnt, tx.gcct.prop, tx.gcct.propv2, gcct.infect.dur,
                 gcct.asympt.tests.py, gcct.asympt.tests,  gccttxpy,
                 syph.incid, syph.pia, syph.nnt, tx.syph.prop, tx.syph.propv2, syph.infect.dur, syph.asympt.tests.py,
                 syph.asympt.tests, syphtxpy, syphearlytxpy, syphlatetxpy,

                 # Group 1
                 gcct.incid.g1, gcct.pia.g1, gcct.nnt.g1,
                 gcct.asympt.tests.py.g1,
                 gcct.asympt.tests.g1, gccttxpy.g1,
                 syph.incid.g1, syph.pia.g1, syph.nnt.g1,
                 syph.asympt.tests.py.g1, syph.asympt.tests.g1,
                 syphtxpy.g1, syphearlytxpy.g1, syphlatetxpy.g1,

                 # Group 2
                 gcct.incid.g2, gcct.pia.g2, gcct.nnt.g2,
                 gcct.asympt.tests.py.g2, gcct.asympt.tests.g2, gccttxpy.g2,
                 syph.incid.g2, syph.pia.g2, syph.nnt.g2, syph.asympt.tests.py.g2, syph.asympt.tests.g2,
                 syphtxpy.g2, syphearlytxpy.g2, syphlatetxpy.g2
)


for (i in seq_along(sims)) {

  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int

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

  # Proportion treated (over last year)
  vec.tx.gcct.prop <- unname(colMeans(tail(sim$epi$tx.gcct.prop, 52)))
  vec.tx.syph.prop <- unname(colMeans(tail(sim$epi$tx.syph.prop, 52)))
  df$tx.gcct.prop[i] <- paste0(round(quantile(vec.tx.gcct.prop, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.tx.gcct.prop, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.tx.gcct.prop, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
  df$tx.syph.prop[i] <- paste0(round(quantile(vec.tx.syph.prop, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.tx.syph.prop, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.tx.syph.prop, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")

  vec.tx.gcct.propv2 <- (unname(colSums(tail(sim$epi$txGC, 52))) + unname(colSums(tail(sim$epi$txCT, 52)))) / (unname(colSums(tail(sim$epi$incid.gc, 52))) + unname(colSums(tail(sim$epi$incid.ct, 52))))
  vec.tx.syph.propv2 <- (unname(colSums(tail(sim$epi$txsyph, 52)))) / unname(colSums(tail(sim$epi$incid.syph, 52)))
  df$tx.gcct.propv2[i] <- paste0(round(quantile(vec.tx.gcct.propv2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                 " (", round(quantile(vec.tx.gcct.propv2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                 " - ", round(quantile(vec.tx.gcct.propv2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                 ")")
  df$tx.syph.propv2[i] <- paste0(round(quantile(vec.tx.syph.propv2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                 " (", round(quantile(vec.tx.syph.propv2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                 " - ", round(quantile(vec.tx.syph.propv2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                 ")")
  # Infection duration until recovery
  vec.gcct.infect.dur <- unname(colMeans(tail(sim$epi$gcct.infect.dur, 52), na.rm = TRUE))
  vec.syph.infect.dur <- unname(colMeans(tail(sim$epi$syph.infect.dur, 52), na.rm = TRUE))
  df$gcct.infect.dur[i] <- paste0(round(quantile(vec.gcct.infect.dur, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                  " (", round(quantile(vec.gcct.infect.dur, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                  " - ", round(quantile(vec.gcct.infect.dur, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                  ")")
  df$syph.infect.dur[i] <- paste0(round(quantile(vec.syph.infect.dur, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                  " (", round(quantile(vec.syph.infect.dur, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                  " - ", round(quantile(vec.syph.infect.dur, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                  ")")


  # PIA (Cumulative)
  # incid.gcct <- unname(colSums(sim$epi$incid.gcct))
  incid.gcct <- unname(colSums(sim$epi$incid.gc)) + unname(colSums(sim$epi$incid.ct))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  #incid.gcct.g1 <- unname(colSums(sim$epi$incid.gcct.tttraj1))
  incid.gcct.g1 <- unname(colSums(sim$epi$incid.gc.tttraj1)) + unname(colSums(sim$epi$incid.ct.tttraj1))
  vec.nia.gcct.g1 <- round(incid.base.gcct.g1 - incid.gcct.g1, 1)
  vec.pia.gcct.g1 <- vec.nia.gcct.g1/incid.gcct.g1

  #incid.gcct.g2 <- unname(colSums(sim$epi$incid.gcct.tttraj2))
  incid.gcct.g2 <- unname(colSums(sim$epi$incid.gcct.tttraj2)) + unname(colSums(sim$epi$incid.ct.tttraj2))
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
  vec.pia.syph.g2 <- vec.nia.syph.g2/incid.syph.g2

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

  gcct.asympt.tests.py <-  ifelse(py > 0, 52 * (gcct.asympt.tests / py), 0)
  gcct.asympt.tests.py.g1 <-  ifelse(py.g1 > 0,52 * (gcct.asympt.tests.g1 / py.g1), 0)
  gcct.asympt.tests.py.g2 <-  ifelse(py > 0,52 * (gcct.asympt.tests.g2 / py.g2), 0)
  syph.asympt.tests.py <-  ifelse(py > 0,52 * (syph.asympt.tests / py), 0)
  syph.asympt.tests.py.g1 <-  ifelse(py.g1 > 0, 52 * (syph.asympt.tests.g1 / py.g1), 0)
  syph.asympt.tests.py.g2 <-  ifelse(py.g2 > 0, 52 * (syph.asympt.tests.g2 / py.g2), 0)

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
  vec.gcct.nnt <- (gcct.asympt.tests - tests.gcct.base) / (incid.base.gcct - incid.gcct)
  vec.gcct.nnt.g1 <- (gcct.asympt.tests.g1 - tests.gcct.base.g1) / (incid.base.gcct.g1 - incid.gcct.g1)
  vec.gcct.nnt.g2 <- (gcct.asympt.tests.g2 - tests.gcct.base.g2) / (incid.base.gcct.g2 - incid.gcct.g2)

  vec.syph.nnt <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - incid.syph)
  vec.syph.nnt.g1 <- (syph.asympt.tests.g1  - tests.syph.base.g1) / (incid.base.syph.g1 - incid.syph.g1)
  vec.syph.nnt.g2 <- (syph.asympt.tests.g2  - tests.syph.base.g2) / (incid.base.syph.g2 - incid.syph.g2)

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

  # Tx per person-year infected
  vec.tx.gcctpy <- unname(colMeans(52 * (sim$epi$txGC + sim$epi$txCT) / (sim$epi$num * (sim$epi$prev.gc + sim$epi$prev.ct))))
  vec.tx.gcctpy.g1 <- unname(colMeans(52 * (sim$epi$txGC.tttraj1 + sim$epi$txCT.tttraj1) / (sim$epi$tt.traj.sti1 * (sim$epi$prev.gc.tttraj1 + sim$epi$prev.ct.tttraj1))))
  vec.tx.gcctpy.g2 <- unname(colMeans(52 * (sim$epi$txGC.tttraj2 + sim$epi$txCT.tttraj2) / (sim$epi$tt.traj.sti2 * (sim$epi$prev.gc.tttraj2 + sim$epi$prev.ct.tttraj2))))
  # vec.tx.gcctpy <- unname(colMeans(52 * (sim$epi$txGC_asympt + sim$epi$txCT_asympt) / (sim$epi$num * (sim$epi$prev.gcct))))
  # vec.tx.gcctpy.g1 <- unname(colMeans(52 * (sim$epi$txGC_asympt.tttraj1 + sim$epi$txCT_asympt.tttraj1) / (sim$epi$num * (sim$epi$prev.gcct.g1))))
  # vec.tx.gcctpy.g2 <- unname(colMeans(52 * (sim$epi$txGC_asympt.tttraj2 + sim$epi$txCT_asympt.tttraj2) / (sim$epi$num * (sim$epi$prev.gcct.g2))))

  df$gccttxpy[i] <- paste0(round(quantile(vec.tx.gcctpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.tx.gcctpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.tx.gcctpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$gccttxpy.g1[i] <- paste0(round(quantile(vec.tx.gcctpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.gcctpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.gcctpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$gccttxpy.g2[i] <- paste0(round(quantile(vec.tx.gcctpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.gcctpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.gcctpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  vec.tx.syphpy <- unname(colMeans(52 * sim$epi$txsyph / (sim$epi$num * sim$epi$prev.syph)))
  vec.tx.syphpy.g1 <- unname(colMeans(52 * sim$epi$txsyph.tttraj1 / (sim$epi$tt.traj.sti1 * sim$epi$prev.syph.tttraj1)))
  vec.tx.syphpy.g2 <- unname(colMeans(52 * sim$epi$txsyph.tttraj2 / (sim$epi$tt.traj.sti2 * sim$epi$prev.syph.tttraj2)))

  # Prev.primsecosyph (Incub, prim, syph)
  df.prev.earlysyph <- sim$epi$txearlysyph[1:521,] / (sim$epi$num[1:521,] * sim$epi$prev.primsecosyph[1:521,])
  df.prev.earlysyph.g1 <- sim$epi$txearlysyph[1:521,] / (sim$epi$tt.traj.sti1[1:521,] * sim$epi$prev.primsecosyph[1:521,])
  df.prev.earlysyph.g2 <- sim$epi$txearlysyph[1:521,] / (sim$epi$tt.traj.sti2[1:521,] * sim$epi$prev.primsecosyph[1:521,])
  # Remove NaNs
  for (j in 1:ncol(df.prev.earlysyph)) {
    df.prev.earlysyph[which(is.nan(df.prev.earlysyph[, j])), j] <- 0.0
  }
  for (j in 1:ncol(df.prev.earlysyph.g1)) {
    df.prev.earlysyph.g1[which(is.nan(df.prev.earlysyph.g1[, j])), j] <- 0.0
  }
  for (j in 1:ncol(df.prev.earlysyph.g2)) {
    df.prev.earlysyph.g2[which(is.nan(df.prev.earlysyph.g2[, j])), j] <- 0.0
  }
  vec.tx.earlysyphpy <- unname(colMeans(52 * df.prev.earlysyph))
  vec.tx.earlysyphpy.g1 <- unname(colMeans(52 * df.prev.earlysyph.g1))
  vec.tx.earlysyphpy.g2 <- unname(colMeans(52 * df.prev.earlysyph.g2))

  # Late syphilis prevalence = Prev.syph - Prev.primsecosyph
  df.prev.latesyph <- sim$epi$txlatesyph[1:521,] / (sim$epi$num[1:521,] * (sim$epi$prev.syph[1:521, ] - sim$epi$prev.primsecosyph[1:521,]))
  df.prev.latesyph.g1 <- sim$epi$txlatesyph[1:521,] / (sim$epi$tt.traj.sti1[1:521,] * (sim$epi$prev.syph[1:521, ] - sim$epi$prev.primsecosyph[1:521,]))
  df.prev.latesyph.g2 <- sim$epi$txlatesyph[1:521,] / (sim$epi$tt.traj.sti2[1:521,] * (sim$epi$prev.syph[1:521, ] - sim$epi$prev.primsecosyph[1:521,]))

  # Remove NaNs
  for (j in 1:ncol(df.prev.latesyph)) {
    df.prev.latesyph[which(is.nan(df.prev.latesyph[, j])), j] <- 0.0
  }
  for (j in 1:ncol(df.prev.latesyph.g1)) {
    df.prev.latesyph.g1[which(is.nan(df.prev.latesyph.g1[, j])), j] <- 0.0
  }
  for (j in 1:ncol(df.prev.latesyph.g2)) {
    df.prev.latesyph.g2[which(is.nan(df.prev.latesyph.g2[, j])), j] <- 0.0
  }
  vec.tx.latesyphpy <- unname(colMeans(52 * df.prev.latesyph))
  vec.tx.latesyphpy.g1 <- unname(colMeans(52 * df.prev.latesyph.g1))
  vec.tx.latesyphpy.g2 <- unname(colMeans(52 * df.prev.latesyph.g2))

  df$syphtxpy[i] <- paste0(round(quantile(vec.tx.syphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.tx.syphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.tx.syphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$syphtxpy.g1[i] <- paste0(round(quantile(vec.tx.syphpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.syphpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.syphpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")
  df$syphtxpy.g2[i] <- paste0(round(quantile(vec.tx.syphpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.tx.syphpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.tx.syphpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  df$syphearlytxpy[i] <- paste0(round(quantile(vec.tx.earlysyphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                " (", round(quantile(vec.tx.earlysyphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                " - ", round(quantile(vec.tx.earlysyphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                ")")
  df$syphearlytxpy.g1[i] <- paste0(round(quantile(vec.tx.earlysyphpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                   " (", round(quantile(vec.tx.earlysyphpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                   " - ", round(quantile(vec.tx.earlysyphpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                   ")")
  df$syphearlytxpy.g2[i] <- paste0(round(quantile(vec.tx.earlysyphpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                   " (", round(quantile(vec.tx.earlysyphpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                   " - ", round(quantile(vec.tx.earlysyphpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                   ")")
  df$syphlatetxpy[i] <- paste0(round(quantile(vec.tx.latesyphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                               " (", round(quantile(vec.tx.latesyphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                               " - ", round(quantile(vec.tx.latesyphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                               ")")
  df$syphlatetxpy.g1[i] <- paste0(round(quantile(vec.tx.latesyphpy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                  " (", round(quantile(vec.tx.latesyphpy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                  " - ", round(quantile(vec.tx.latesyphpy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                  ")")
  df$syphlatetxpy.g2[i] <- paste0(round(quantile(vec.tx.latesyphpy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                  " (", round(quantile(vec.tx.latesyphpy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                  " - ", round(quantile(vec.tx.latesyphpy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                  ")")

  cat("*")

}

df

write.csv(df, "analysis/New STD Table 2.csv")
