## STI Testing Guidelines Table 3
# Sensitivity Analyses in appendix - STI specific

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
# Reference scenario here
load("data/followup/sim.n3191.rda")
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

## -
# Baseline sim = 0.2 coverage of annual and HR

# Screening Intervals:
# 3191 as reference
# 3189, 3190, 3191 (ref), 3192, 3193: Annual = 182 days, 273 days, 364 days (ref), 448 days, 539 days, HR = 50%, ANN = 50%, 182 days
# 3194-3195, 3197-3198: Higher-risk = 28 days, 91 days, 182 days (ref), 273 days, 364 days, HR = 50%, Ann = 50%, 364 days
#
# Treatment Completion
# 3199, 3204, 3209, 3214: Annual = 0.0 - 1.0 (ref) by 0.25, 364 days, HR = 0%, 182 days
#
# Partner Cutoff for Higher-Risk
# 3221:3229 Higher-risk = 1 (ref) to 10 by 1

# Newer way:
sims <- c(3189:3190, 3191, 3192:3193, 3194:3195, 3191, 3197:3198, 3199, 3204, 3209, 3214, 3191, 3191, 3221:3229)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
partcut <- rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
gc.tests <- rep(NA, length(sims))
gctx <- rep(NA, length(sims))
gctxpy <- rep(NA, length(sims))

gc.incid.g1 <- rep(NA, length(sims))
gc.pia.g1 <- rep(NA, length(sims))
gc.tests.g1 <- rep(NA, length(sims))
gctx.g1 <- rep(NA, length(sims))
gctxpy.g1 <- rep(NA, length(sims))

gc.incid.g2 <- rep(NA, length(sims))
gc.pia.g2 <- rep(NA, length(sims))
gc.tests.g2 <- rep(NA, length(sims))
gctx.g2 <- rep(NA, length(sims))
gctxpy.g2 <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
ct.tests <- rep(NA, length(sims))
cttx <- rep(NA, length(sims))
cttxpy <- rep(NA, length(sims))

ct.incid.g1 <- rep(NA, length(sims))
ct.pia.g1 <- rep(NA, length(sims))
ct.tests.g1 <- rep(NA, length(sims))
cttx.g1 <- rep(NA, length(sims))
cttxpy.g1 <- rep(NA, length(sims))

ct.incid.g2 <- rep(NA, length(sims))
ct.pia.g2 <- rep(NA, length(sims))
ct.tests.g2 <- rep(NA, length(sims))
cttx.g2 <- rep(NA, length(sims))
cttxpy.g2 <- rep(NA, length(sims))

syph.incid <- rep(NA, length(sims))
syph.pia <- rep(NA, length(sims))
syph.tests <- rep(NA, length(sims))
syphearlytx2 <- rep(NA, length(sims))
syphlatetx <- rep(NA, length(sims))
syphtxpy <- rep(NA, length(sims))

syph.incid.g1 <- rep(NA, length(sims))
syph.pia.g1 <- rep(NA, length(sims))
syph.tests.g1 <- rep(NA, length(sims))
syphearlytx.g1 <- rep(NA, length(sims))
syphlatetx.g1 <- rep(NA, length(sims))
syphtxpy.g1 <- rep(NA, length(sims))

syph.incid.g2 <- rep(NA, length(sims))
syph.pia.g2 <- rep(NA, length(sims))
syph.tests.g2 <- rep(NA, length(sims))
syphtx.g2 <- rep(NA, length(sims))
syphearlytx.g2 <- rep(NA, length(sims))
syphlatetx.g2 <- rep(NA, length(sims))
syphtxpy.g2 <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, partcut,
                 # Overall
                 gc.incid, gc.tests, gc.pia, gctx, gctxpy,
                 ct.incid, ct.tests, ct.pia, cttx, cttxpy,
                 syph.incid, syph.tests, syph.pia, syphtx, earlysyphtxpy, latesyphtxpy,

                 # Group 1
                 gc.incid.g1, gc.tests.g1, gc.pia.g1, gctx.g1, gctxpy.g1,
                 ct.incid.g1, ct.tests.g1, ct.pia.g1, cttx.g1, cttxpy.g1,
                 syph.incid.g1, syph.tests.g1, syph.pia.g1, syphtx.g1, earlysyphtxpy.g1, latesyphtxpy.g1,

                 # Group 2
                 gc.incid.g2, gc.tests.g2, gc.pia.g2, gctx.g2, gctxpy.g2,
                 ct.incid.g2, ct.tests.g2, ct.pia.g2, cttx.g2, cttxpy.g2,
                 syph.incid.g2, syph.tests.g2, syph.pia.g2, syphtx.g2, earlysyphtxpy.g2, latesyphtxpy.g2)


for (i in seq_along(sims)) {

  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int
  df$partcut[i] <- sim$param$partnercut

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
  ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc, na.rm = TRUE)) * 1000
  vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
  vec.pia.gc <- vec.nia.gc/ir.base.gc
  vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]

  ir.comp.gc.g1 <- unname(colMeans(sim$epi$ir100.gc.tttraj1, na.rm = TRUE)) * 1000
  vec.nia.gc.g1 <- round(ir.base.gc.g1 - ir.comp.gc.g1, 1)
  vec.pia.gc.g1 <- vec.nia.gc.g1/ir.base.gc
  vec.pia.gc.g1 <- vec.pia.gc.g1[vec.pia.gc.g1 > -Inf]

  ir.comp.gc.g2 <- unname(colMeans(sim$epi$ir100.gc.tttraj2, na.rm = TRUE)) * 1000
  vec.nia.gc.g2 <- round(ir.base.gc.g2 - ir.comp.gc.g2, 1)
  vec.pia.gc.g2 <- vec.nia.gc.g2/ir.base.gc
  vec.pia.gc.g2 <- vec.pia.gc.g2[vec.pia.gc.g2 > -Inf]

  ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct, na.rm = TRUE)) * 1000
  vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
  vec.pia.ct <- vec.nia.ct/ir.base.ct
  vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]

  ir.comp.ct.g1 <- unname(colMeans(sim$epi$ir100.ct.tttraj1, na.rm = TRUE)) * 1000
  vec.nia.ct.g1 <- round(ir.base.ct.g1 - ir.comp.ct.g1, 1)
  vec.pia.ct.g1 <- vec.nia.ct.g1/ir.base.ct
  vec.pia.ct.g1 <- vec.pia.ct.g1[vec.pia.ct.g1 > -Inf]

  ir.comp.ct.g2 <- unname(colMeans(sim$epi$ir100.ct.tttraj2, na.rm = TRUE)) * 1000
  vec.nia.ct.g2 <- round(ir.base.ct.g2 - ir.comp.ct.g2, 1)
  vec.pia.ct.g2 <- vec.nia.ct.g2/ir.base.ct
  vec.pia.ct.g2 <- vec.pia.ct.g2[vec.pia.ct.g2 > -Inf]

  ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph, na.rm = TRUE)) * 1000
  vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
  vec.pia.syph <- vec.nia.syph/ir.base.syph
  vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]

  ir.comp.syph.g1 <- unname(colMeans(sim$epi$ir100.syph.tttraj1, na.rm = TRUE)) * 1000
  vec.nia.syph.g1 <- round(ir.base.syph.g1 - ir.comp.syph.g1, 1)
  vec.pia.syph.g1 <- vec.nia.syph.g1/ir.base.syph
  vec.pia.syph.g1 <- vec.pia.syph.g1[vec.pia.syph.g1 > -Inf]

  ir.comp.syph.g2 <- unname(colMeans(sim$epi$ir100.syph.tttraj2, na.rm = TRUE)) * 1000
  vec.nia.syph.g2 <- round(ir.base.syph.g2 - ir.comp.syph.g2, 1)
  vec.pia.syph.g2 <- vec.nia.syph.g2/ir.base.syph
  vec.pia.syph.g2 <- vec.pia.syph.g2[vec.pia.syph.g2 > -Inf]

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

  #Tests (over first year or cumulative?)
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

  df$gc.tests[i] <- paste0(round(quantile(gc.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(gc.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(gc.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
  df$gc.tests.g1[i] <- paste0(round(quantile(gc.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(gc.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(gc.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")
  df$gc.tests.g2[i] <- paste0(round(quantile(gc.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                               " (", round(quantile(gc.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                               " - ", round(quantile(gc.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                               ")")

  df$ct.tests[i] <- paste0(round(quantile(ct.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(ct.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(ct.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
  df$ct.tests.g1[i] <- paste0(round(quantile(ct.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                              " (", round(quantile(ct.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                              " - ", round(quantile(ct.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                              ")")
  df$ct.tests.g2[i] <- paste0(round(quantile(ct.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                              " (", round(quantile(ct.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                              " - ", round(quantile(ct.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                              ")")

  df$syph.tests[i] <- paste0(round(quantile(syph.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                           " (", round(quantile(syph.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                           " - ", round(quantile(syph.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                           ")")
  df$syph.tests.g1[i] <- paste0(round(quantile(syph.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                              " (", round(quantile(syph.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                              " - ", round(quantile(syph.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                              ")")
  df$syph.tests.g2[i] <- paste0(round(quantile(syph.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                              " (", round(quantile(syph.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                              " - ", round(quantile(syph.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                              ")")

  # Treated in first year or Cumulative?
  # vec.tx.gc <- unname(colSums(sim$epi$txGC))
  # vec.tx.gc.g1 <- unname(colSums(sim$epi$txGC.tttraj1))
  # vec.tx.gc.g2 <- unname(colSums(sim$epi$txGC.tttraj2))
  vec.tx.gc <- as.numeric(colSums(head(sim$epi$txGC, 52)))
  vec.tx.gc.g1 <- as.numeric(colSums(head(sim$epi$txGC.tttraj1, 52)))
  vec.tx.gc.g2 <- as.numeric(colSums(head(sim$epi$txGC.tttraj2, 52)))

  # vec.tx.ct <- unname(colSums(sim$epi$txCT))
  # vec.tx.ct.g1 <- unname(colSums(sim$epi$txCT.tttraj1))
  # vec.tx.ct.g2 <- unname(colSums(sim$epi$txCT.tttraj2))
  vec.tx.ct <- as.numeric(colSums(head(sim$epi$txCT, 52)))
  vec.tx.ct.g1 <- as.numeric(colSums(head(sim$epi$txCT.tttraj1, 52)))
  vec.tx.ct.g2 <- as.numeric(colSums(head(sim$epi$txCT.tttraj2, 52)))

  # vec.tx.syph <- unname(colSums(sim$epi$txsyph))
  # vec.tx.syph.g1 <- unname(colSums(sim$epi$txsyph.tttraj1))
  # vec.tx.syph.g2 <- unname(colSums(sim$epi$txsyph.tttraj2))
  vec.tx.syph <- as.numeric(colSums(head(sim$epi$txsyph, 52)))
  vec.tx.syph.g1 <- as.numeric(colSums(head(sim$epi$txsyph.tttraj1, 52)))
  vec.tx.syph.g2 <- as.numeric(colSums(head(sim$epi$txsyph.tttraj2, 52)))

  # vec.tx.earlysyph <- unname(colSums(sim$epi$txearlysyph))
  # vec.tx.earlysyph.g1 <- unname(colSums(sim$epi$txearlysyph.tttraj1))
  # vec.tx.earlysyph.g2 <- unname(colSums(sim$epi$txearlysyph.tttraj2))
  vec.tx.earlysyph <- as.numeric(colSums(head(sim$epi$txearlysyph, 52)))
  vec.tx.earlysyph.g1 <- as.numeric(colSums(head(sim$epi$txearlysyph.tttraj1, 52)))
  vec.tx.earlysyph.g2 <- as.numeric(colSums(head(sim$epi$txearlysyph.tttraj2, 52)))

  # vec.tx.latesyph <- unname(colSums(sim$epi$txlatesyph))
  # vec.tx.latesyph.g1 <- unname(colSums(sim$epi$txlatesyph.tttraj1))
  # vec.tx.latesyph.g2 <- unname(colSums(sim$epi$txlatesyph.tttraj2))
  vec.tx.latesyph <- as.numeric(colSums(head(sim$epi$txlatesyph, 52)))
  vec.tx.latesyph.g1 <- as.numeric(colSums(head(sim$epi$txlatesyph.tttraj1, 52)))
  vec.tx.latesyph.g2 <- as.numeric(colSums(head(sim$epi$txlatesyph.tttraj2, 52)))

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

  df$syphtx[i] <- paste0(round(quantile(vec.tx.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                       " (", round(quantile(vec.tx.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                       " - ", round(quantile(vec.tx.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                       ")")
  df$syphtx.g1[i] <- paste0(round(quantile(vec.tx.syph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                          " (", round(quantile(vec.tx.syph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                          " - ", round(quantile(vec.tx.syph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                          ")")
  df$syphtx.g2[i] <- paste0(round(quantile(vec.tx.syph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                          " (", round(quantile(vec.tx.syph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                          " - ", round(quantile(vec.tx.syph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                          ")")

  df$earlysyphtx[i] <- paste0(round(quantile(vec.tx.earlysyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.tx.earlysyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.tx.earlysyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
  df$earlysyphtx.g1[i] <- paste0(round(quantile(vec.tx.earlysyph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.tx.earlysyph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.tx.earlysyph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
  df$earlysyphtx.g2[i] <- paste0(round(quantile(vec.tx.earlysyph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.tx.earlysyph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.tx.earlysyph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")

  df$latesyphtx[i] <- paste0(round(quantile(vec.tx.latesyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                         " (", round(quantile(vec.tx.latesyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                         " - ", round(quantile(vec.tx.latesyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                         ")")
  df$latesyphtx.g1[i] <- paste0(round(quantile(vec.tx.latesyph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.tx.latesyph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.tx.latesyph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")
  df$latesyphtx.g2[i] <- paste0(round(quantile(vec.tx.latesyph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
                            " (", round(quantile(vec.tx.latesyph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
                            " - ", round(quantile(vec.tx.latesyph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
                            ")")


  #txGC or txGC_asympt? (Over first year or cumulative?)
  #Update n values to num.tttraj1 and num.tttraj2
  # treated over first year?
  vec.tx.gcpy <- unname(colMeans(52 * sim$epi$txGC / (sim$epi$num * sim$epi$prev.gc)))
  vec.tx.gcpy.g1 <- unname(colMeans(52 * sim$epi$txGC.tttraj1 / (sim$epi$tt.traj.gc1 * sim$epi$prev.gc.tttraj1)))
  vec.tx.gcpy.g2 <- ifelse(median(unname(colMeans(sim$epi$tt.gc.sti2))) > 0,
                           unname(colMeans(52 * sim$epi$txGC.tttraj2 / (sim$epi$tt.traj.gc2 * sim$epi$prev.gc.tttraj2))),
                           0)
  vec.tx.ctpy <- unname(colMeans(52 * sim$epi$txCT / (sim$epi$num * sim$epi$prev.ct)))
  vec.tx.ctpy.g1 <- unname(colMeans(52 * sim$epi$txCT.tttraj1 / (sim$epi$tt.traj.ct1 * sim$epi$prev.ct.tttraj1)))
  vec.tx.ctpy.g2 <- ifelse(median(unname(colMeans(sim$epi$tt.ct.sti2))) > 0,
                           unname(colMeans(52 * sim$epi$txCT.tttraj2 / (sim$epi$tt.traj.ct2 * sim$epi$prev.ct.tttraj2))),
                           0)
  vec.tx.syphpy <- unname(colMeans(52 * sim$epi$txsyph / (sim$epi$num * sim$epi$prev.syph)))
  vec.tx.syphpy.g1 <- unname(colMeans(52 * sim$epi$txsyph.tttraj1 / (sim$epi$tt.traj.syph1 * sim$epi$prev.syph.tttraj1)))
  vec.tx.syphpy.g2 <- ifelse(median(unname(colMeans(sim$epi$tt.syph.sti2))) > 0,
                           unname(colMeans(52 * sim$epi$txsyph.tttraj2 / (sim$epi$tt.traj.syph2 * sim$epi$prev.syph.tttraj2))),
                           0)


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


  cat("*")

}

df

write.csv(df, "analysis/STD Table 3.csv")
