## STI Testing Guidelines Table 2
# Sensitivity Analyses

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
# Reference scenario here
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
ir.base.sti.g1 <- unname(colMeans(sim.base$epi$ir100.sti.tttraj1)) * 1000
ir.base.sti.g2 <- unname(colMeans(sim.base$epi$ir100.sti.tttraj2)) * 1000


## -
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
sims <- c(3000, 3189:3193, 3194:3198, 3221:3229)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
partcut <- rep(NA, length(sims))

sti.incid <- rep(NA, length(sims))
sti.pia <- rep(NA, length(sims))
sti.tests <- rep(NA, length(sims))
tx <- rep(NA, length(sims))
txperpy <- rep(NA, length(sims))

sti.incid.g1 <- rep(NA, length(sims))
sti.pia.g1 <- rep(NA, length(sims))
sti.tests.g1 <- rep(NA, length(sims))
tx.g1 <- rep(NA, length(sims))
txperpy.g1 <- rep(NA, length(sims))

sti.incid.g2 <- rep(NA, length(sims))
sti.pia.g2 <- rep(NA, length(sims))
sti.tests.g2 <- rep(NA, length(sims))
tx.g2 <- rep(NA, length(sims))
txperpy.g2 <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, partcut,
                 # Overall
                 sti.incid,  sti.pia, sti.tests,tx, txperpy,
                 # Group 1
                 sti.incid.g1, sti.pia.g1, sti.tests.g1,  tx.g1, txperpy.g1,
                 # Group 2
                 sti.incid.g2, sti.pia.g2, sti.tests.g2,  tx.g2, txperpy.g2
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
  vec.ir.sti <- unname(colMeans(tail(sim$epi$ir100.sti, 52)))
  vec.ir.sti.g1 <- unname(colMeans(tail(sim$epi$ir100.sti.tttraj1, 52)))
  vec.ir.sti.g2 <- unname(colMeans(tail(sim$epi$ir100.sti.tttraj2, 52)))
  df$sti.incid[i] <- paste0(round(quantile(vec.ir.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$sti.incid.g1[i] <- paste0(round(quantile(vec.ir.sti.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.sti.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.sti.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$sti.incid.g2[i] <- paste0(round(quantile(vec.ir.sti.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.ir.sti.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.ir.sti.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  # PIA (Cumulative)
  ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti, na.rm = TRUE)) * 1000
  vec.nia.sti <- round(ir.base.sti - ir.comp.sti, 1)
  vec.pia.sti <- vec.nia.sti/ir.base.sti
  vec.pia.sti <- vec.pia.sti[vec.pia.sti > -Inf]

  ir.comp.sti.g1 <- unname(colMeans(sim$epi$ir100.sti.tttraj1, na.rm = TRUE)) * 1000
  vec.nia.sti.g1 <- round(ir.base.sti.g1 - ir.comp.sti.g1, 1)
  vec.pia.sti.g1 <- vec.nia.sti.g1/ir.base.sti
  vec.pia.sti.g1 <- vec.pia.sti.g1[vec.pia.sti.g1 > -Inf]

  ir.comp.sti.g2 <- unname(colMeans(sim$epi$ir100.sti.tttraj2, na.rm = TRUE)) * 1000
  vec.nia.sti.g2 <- round(ir.base.sti.g2 - ir.comp.sti.g2, 1)
  vec.pia.sti.g2 <- vec.nia.sti.g2/ir.base.sti
  vec.pia.sti.g2 <- vec.pia.sti.g2[vec.pia.sti.g2 > -Inf]

  df$sti.pia[i] <- paste0(round(quantile(vec.pia.sti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.sti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.sti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")
  df$sti.pia.g1[i] <- paste0(round(quantile(vec.pia.sti.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.sti.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.sti.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$sti.pia.g2[i] <- paste0(round(quantile(vec.pia.sti.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.sti.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.sti.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")

  #Tests (Cumulative or over first year?)
  sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))
  sti.sympt.tests <- unname(colSums(sim$epi$stisympttests, na.rm = TRUE))
  sti.tests <- sti.asympt.tests + sti.sympt.tests

  sti.asympt.tests.g1 <- unname(colSums(sim$epi$stiasympttests.tttraj1, na.rm = TRUE))
  sti.sympt.tests.g1 <- unname(colSums(sim$epi$stisympttests.tttraj1, na.rm = TRUE))
  sti.tests.g1 <- sti.asympt.tests.g1 + sti.sympt.tests.g1

  sti.asympt.tests.g2 <- unname(colSums(sim$epi$stiasympttests.tttraj2, na.rm = TRUE))
  sti.sympt.tests.g2 <- unname(colSums(sim$epi$stisympttests.tttraj2, na.rm = TRUE))
  sti.tests.g2 <- sti.asympt.tests.g2 + sti.sympt.tests.g2

  df$sti.tests[i] <- paste0(round(quantile(sti.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(sti.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(sti.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
  df$sti.tests.g1[i] <- paste0(round(quantile(sti.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(sti.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(sti.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
  df$sti.tests.g2[i] <- paste0(round(quantile(sti.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(sti.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(sti.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")


  vec.tx <- unname(colSums(sim$epi$txSTI))
  vec.tx.g1 <- unname(colSums(sim$epi$txSTI.tttraj1))
  vec.tx.g2 <- unname(colSums(sim$epi$txSTI.tttraj2))
  df$tx[i] <- paste0(round(quantile(vec.tx, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(vec.tx, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(vec.tx, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
  df$tx.g1[i] <- paste0(round(quantile(vec.tx.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                     " (", round(quantile(vec.tx.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                     " - ", round(quantile(vec.tx.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                     ")")
  df$tx.g2[i] <- paste0(round(quantile(vec.tx.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                     " (", round(quantile(vec.tx.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                     " - ", round(quantile(vec.tx.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                     ")")

  #txSTI or txSTI_asympt? (Over first year or cumulative?)
  #Update n values to num.tttraj1 and num.tttraj2
  vec.tx.stipy <- unname(colMeans(52 * sim$epi$txSTI / (sim$epi$num * sim$epi$prev.sti)))
  vec.tx.stipy.g1 <- unname(colMeans(52 * sim$epi$txSTI.tttraj1 / (sim$epi$tt.traj.sti1 * sim$epi$prev.sti.tttraj1)))
  df.prev.sti.tttraj2 <- sim$epi$txSTI.tttraj2[1:521,] / (sim$epi$tt.traj.sti2[1:521,] * sim$epi$prev.sti.tttraj2[1:521,])

  # Remove NaNs
  for (j in 1:ncol(df.prev.sti.tttraj2)) {

        df.prev.sti.tttraj2[which(is.nan(df.prev.sti.tttraj2[, j])), j] <- 0.0

  }
  vec.tx.stipy.g2 <- unname(colMeans(52 * df.prev.sti.tttraj2))
  df$txperpy[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

  df$txperpy.g1[i] <- paste0(round(quantile(vec.tx.stipy.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.stipy.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.stipy.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  df$txperpy.g2[i] <- paste0(round(quantile(vec.tx.stipy.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                            " (", round(quantile(vec.tx.stipy.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                            " - ", round(quantile(vec.tx.stipy.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                            ")")

  cat("*")

}

df

write.csv(df, "analysis/STD Table 2.csv")

