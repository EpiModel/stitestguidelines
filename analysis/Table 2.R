## STI Testing Guidelines Table 2
# Sensitivity Analyses

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
# Reference scenario here
load("data/followup/sim.n3000.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

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

## -
# Baseline compared to 5% HR coverage

# Screening Intervals:
# 3191 as reference
# 3189, 3190, 3191 (ref), 3192, 3193: Annual = 182 days, 273 days, 364 days (ref), 448 days, 539 days, HR = 50%, ANN = 50%, 182 days
# 3194-3195, 3197-3198: Higher-risk = 28 days, 91 days, 182 days (ref), 273 days, 364 days, HR = 50%, Ann = 50%, 364 days
#
# Treatment Completion - for appendix
# 3199, 3204, 3209, 3214: Annual = 0.0 - 1.0 (ref) by 0.25, 364 days, HR = 0%, 182 days
#
# Partner Cutoff for Higher-Risk
# 3221:3229 Higher-risk = 1 (ref) to 10 by 1

# Newer way:
sims <- c(3189:3190, 3191, 3192:3193, 3194:3195, 3191, 3197:3198, 3191, 3221:3229)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
stiasympttx <- rep(NA, length(sims))
partcut <- rep(NA, length(sims))

# hiv.incid <- rep(NA, length(sims))
# hiv.incid.low <- rep(NA, length(sims))
# hiv.incid.high <- rep(NA, length(sims))
#
# hiv.hr <- rep(NA, length(sims))
# hiv.hr.low <- rep(NA, length(sims))
# hiv.hr.high <- rep(NA, length(sims))
#
# hiv.pia <- rep(NA, length(sims))
# hiv.pia.low <- rep(NA, length(sims))
# hiv.pia.high <- rep(NA, length(sims))
#
# hiv.nnt <- rep(NA, length(sims))
# hiv.nnt.low <- rep(NA, length(sims))
# hiv.nnt.high <- rep(NA, length(sims))
#
# gc.incid <- rep(NA, length(sims))
# gc.incid.low <- rep(NA, length(sims))
# gc.incid.high <- rep(NA, length(sims))
#
# gc.hr <- rep(NA, length(sims))
# gc.hr.low <- rep(NA, length(sims))
# gc.hr.high <- rep(NA, length(sims))
#
# gc.pia <- rep(NA, length(sims))
# gc.pia.low <- rep(NA, length(sims))
# gc.pia.high <- rep(NA, length(sims))
#
# gc.nnt <- rep(NA, length(sims))
# gc.nnt.low <- rep(NA, length(sims))
# gc.nnt.high <- rep(NA, length(sims))
#
# ct.incid <- rep(NA, length(sims))
# ct.incid.low <- rep(NA, length(sims))
# ct.incid.high <- rep(NA, length(sims))
#
# ct.hr <- rep(NA, length(sims))
# ct.hr.low <- rep(NA, length(sims))
# ct.hr.high <- rep(NA, length(sims))
#
# ct.pia <- rep(NA, length(sims))
# ct.pia.low <- rep(NA, length(sims))
# ct.pia.high <- rep(NA, length(sims))
#
# ct.nnt <- rep(NA, length(sims))
# ct.nnt.low <- rep(NA, length(sims))
# ct.nnt.high <- rep(NA, length(sims))
#
# syph.incid <- rep(NA, length(sims))
# syph.incid.low <- rep(NA, length(sims))
# syph.incid.high <- rep(NA, length(sims))
#
# syph.hr <- rep(NA, length(sims))
# syph.hr.low <- rep(NA, length(sims))
# syph.hr.high <- rep(NA, length(sims))
#
# syph.pia <- rep(NA, length(sims))
# syph.pia.low <- rep(NA, length(sims))
# syph.pia.high <- rep(NA, length(sims))
#
# syph.nnt <- rep(NA, length(sims))
# syph.nnt.low <- rep(NA, length(sims))
# syph.nnt.high <- rep(NA, length(sims))

sti.incid <- rep(NA, length(sims))
sti.incid.low <- rep(NA, length(sims))
sti.incid.high <- rep(NA, length(sims))

sti.hr <- rep(NA, length(sims))
sti.hr.low <- rep(NA, length(sims))
sti.hr.high <- rep(NA, length(sims))

sti.pia <- rep(NA, length(sims))
sti.pia.low <- rep(NA, length(sims))
sti.pia.high <- rep(NA, length(sims))

sti.nnt <- rep(NA, length(sims))
sti.nnt.low <- rep(NA, length(sims))
sti.nnt.high <- rep(NA, length(sims))

# asympt.tx <- rep(NA, length(sims))
# asympt.tx.low <- rep(NA, length(sims))
# asympt.tx.high <- rep(NA, length(sims))
#
# rect.tx <- rep(NA, length(sims))
# rect.tx.low <- rep(NA, length(sims))
# rect.tx.high <- rep(NA, length(sims))
#
# rect.tx <- rep(NA, length(sims))
# rect.tx.low <- rep(NA, length(sims))
# rect.tx.high <- rep(NA, length(sims))

txasympt <- rep(NA, length(sims))
txasympt.low <- rep(NA, length(sims))
txasympt.high <- rep(NA, length(sims))

txCT <- rep(NA, length(sims))
txCT.low <- rep(NA, length(sims))
txCT.high <- rep(NA, length(sims))

txGC <- rep(NA, length(sims))
txGC.low <- rep(NA, length(sims))
txGC.high <- rep(NA, length(sims))

txsyph <- rep(NA, length(sims))
txsyph.low <- rep(NA, length(sims))
txsyph.high <- rep(NA, length(sims))

txearlysyph <- rep(NA, length(sims))
txearlysyph.low <- rep(NA, length(sims))
txearlysyph.high <- rep(NA, length(sims))

txlatesyph <- rep(NA, length(sims))
txlatesyph.low <- rep(NA, length(sims))
txlatesyph.high <- rep(NA, length(sims))

df <- data.frame(anncov, hrcov, annint, hrint, stiasympttx, partcut,
                 # hiv.incid.low, hiv.incid, hiv.incid.high, hiv.hr.low, hiv.hr, hiv.hr.high,
                 # hiv.pia.low, hiv.pia, hiv.pia.high, hiv.nnt.low, hiv.nnt, hiv.nnt.high,
                 # gc.incid.low, gc.incid, gc.incid.high, gc.hr.low, gc.hr, gc.hr.high,
                 # gc.pia.low, gc.pia, gc.pia.high, gc.nnt.low, gc.nnt, gc.nnt.high,
                 # ct.incid.low, ct.incid, ct.incid.high, ct.hr.low, ct.hr, ct.hr.high,
                 # ct.pia.low, ct.pia, ct.pia.high, ct.nnt.low, ct.nnt, ct.nnt.high,
                 # syph.incid.low, syph.incid, syph.incid.high, syph.hr.low, syph.hr, syph.hr.high,
                 # syph.pia.low, syph.pia, syph.pia.high, syph.nnt.low, syph.nnt, syph.nnt.high,

                 # Overall
                 sti.incid.low, sti.incid, sti.incid.high, #sti.hr.low, sti.hr, sti.hr.high,
                 sti.asympttests.low, sti.asympttests, sti.asympttests.high,
                 sti.sympttests.low, sti.sympttests, sti.sympttests.high,
                 sti.pia.low, sti.pia, sti.pia.high, #sti.nnt.low, sti.nnt, sti.nnt.high,
                 # asympt.tx.low, asympt.tx, asympt.tx.high, rect.tx.low, rect.tx, rect.tx.high,
                 txasympt.low, txasympt, txasympt.high, txGC.low, txGC, txGC.high,
                 txCT.low, txCT, txCT.high, txsyph.low, txsyph, txsyph.high,
                 txearlysyph.low, txearlysyph, txearlysyph.high,
                 txlatesyph.low, txlatesyph, txlatesyph.high)

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
  df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
  df$annint[i] <- sim$param$stitest.active.int
  df$hrint[i] <- sim$param$sti.highrisktest.int
  df$stiasympttx[i] <- sim$param$syph.prim.asympt.prob.tx
  df$partcut[i] <- sim$param$partnercut

  # Incidence Rate
  vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
  vec.ir.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
  vec.ir.sti <- unname(colMeans(tail(sim$epi$ir100.sti, 52)))

  # df$hiv.incid.low[i] <- quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$hiv.incid[i] <- quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$hiv.incid.high[i] <- quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$gc.incid.low[i] <- quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$gc.incid[i] <- quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$gc.incid.high[i] <- quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$ct.incid.low[i] <- quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$ct.incid[i] <- quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$ct.incid.high[i] <- quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$syph.incid.low[i] <- quantile(vec.ir.syph, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$syph.incid[i] <- quantile(vec.ir.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$syph.incid.high[i] <- quantile(vec.ir.syph, probs = qnt.high, na.rm = TRUE, names = FALSE)

  df$sti.incid.low[i] <- quantile(vec.ir.sti, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$sti.incid[i] <- quantile(vec.ir.sti, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$sti.incid.high[i] <- quantile(vec.ir.sti, probs = qnt.high, na.rm = TRUE, names = FALSE)

  # # HR
  # num.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
  # denom.hiv <- unname(colMeans(sim.base$epi$ir100)) * 1000
  # hr.vec.hiv <- num.hiv/denom.hiv
  # hr.vec.hiv <- hr.vec.hiv[hr.vec.hiv < Inf]
  # df$hiv.hr.low[i] <- quantile(hr.vec.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$hiv.hr[i] <- quantile(hr.vec.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$hiv.hr.high[i] <- quantile(hr.vec.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # num.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  # denom.gc <- unname(colMeans(tail(sim.base$epi$ir100.gc, 52)))
  # hr.vec.gc <- num.gc/denom.gc
  # hr.vec.gc <- hr.vec.gc[hr.vec.gc < Inf]
  # df$gc.hr.low[i] <- quantile(hr.vec.gc, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$gc.hr[i] <- quantile(hr.vec.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$gc.hr.high[i] <- quantile(hr.vec.gc, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # num.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
  # denom.ct <- unname(colMeans(tail(sim.base$epi$ir100.ct, 52)))
  # hr.vec.ct <- num.ct/denom.ct
  # hr.vec.ct <- hr.vec.ct[hr.vec.ct < Inf]
  # df$ct.hr.low[i] <- quantile(hr.vec.ct, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$ct.hr[i] <- quantile(hr.vec.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$ct.hr.high[i] <- quantile(hr.vec.ct, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # num.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
  # denom.syph <- unname(colMeans(tail(sim.base$epi$ir100.syph, 52)))
  # hr.vec.syph <- num.syph/denom.syph
  # hr.vec.syph <- hr.vec.syph[hr.vec.syph < Inf]
  # df$syph.hr.low[i] <- quantile(hr.vec.syph, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$syph.hr[i] <- quantile(hr.vec.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$syph.hr.high[i] <- quantile(hr.vec.syph, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  #   num.sti <- unname(colMeans(tail(sim$epi$ir100.sti, 52)))
  #   denom.sti <- unname(colMeans(tail(sim.base$epi$ir100.sti, 52)))
  #   hr.vec.sti <- num.sti/denom.sti
  #   hr.vec.sti <- hr.vec.sti[hr.vec.sti < Inf]
  #   df$sti.hr.low[i] <- quantile(hr.vec.sti, probs = qnt.low, na.rm = TRUE, names = FALSE)
  #   df$sti.hr[i] <- quantile(hr.vec.sti, probs = 0.50, na.rm = TRUE, names = FALSE)
  #   df$sti.hr.high[i] <- quantile(hr.vec.sti, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  #PIA
  # ir.comp <- unname(colMeans(sim$epi$ir100)) * 1000
  # vec.nia.hiv <- round(ir.base - ir.comp, 1)
  # vec.pia.hiv <- vec.nia.hiv/ir.base
  # vec.pia.hiv <- vec.pia.hiv[vec.pia.hiv > -Inf]
  # df$hiv.pia.low[i] <- quantile(vec.pia.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$hiv.pia[i] <- quantile(vec.pia.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$hiv.pia.high[i] <- quantile(vec.pia.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  # vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
  # vec.pia.gc <- vec.nia.gc/ir.base.gc
  # vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
  # df$gc.pia.low[i] <- quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$gc.pia[i] <- quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$gc.pia.high[i] <- quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
  # vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
  # vec.pia.ct <- vec.nia.ct/ir.base.ct
  # vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
  # df$ct.pia.low[i] <- quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$ct.pia[i] <- quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$ct.pia.high[i] <- quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
  # vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
  # vec.pia.syph <- vec.nia.syph/ir.base.syph
  # vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
  # df$syph.pia.low[i] <- quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$syph.pia[i] <- quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$syph.pia.high[i] <- quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE)

  ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti)) * 1000
  vec.nia.sti <- round(ir.base.sti - ir.comp.sti, 1)
  vec.pia.sti <- vec.nia.sti/ir.base.sti
  vec.pia.sti <- vec.pia.sti[vec.pia.sti > -Inf]
  df$sti.pia.low[i] <- quantile(vec.pia.sti, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$sti.pia[i] <- quantile(vec.pia.sti, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$sti.pia.high[i] <- quantile(vec.pia.sti, probs = qnt.high, na.rm = TRUE, names = FALSE)

  #NNT
  hiv.tests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
  sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))

  gc.sympt.tests <- unname(colSums(sim$epi$GCsympttests, na.rm = TRUE))
  ct.sympt.tests <- unname(colSums(sim$epi$CTsympttests, na.rm = TRUE))
  syph.sympt.tests <- unname(colSums(sim$epi$syphsympttests, na.rm = TRUE))
  sti.sympt.tests <- unname(colSums(sim$epi$stisympttests, na.rm = TRUE))

  #Tests
  df$sti.asympttests.low[i] <- quantile(sti.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$sti.asympttests[i] <- quantile(sti.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$sti.asympttests.high[i] <- quantile(sti.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE)
  df$sti.sympttests.low[i] <- quantile(sti.sympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$sti.sympttests[i] <- quantile(sti.sympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$sti.sympttests.high[i] <- quantile(sti.sympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE)

  #HIV could be HIV tests or total STI tests
  # vec.hivonly.nnt <- (hiv.tests) / (median(incid.base) - unname(colSums(sim$epi$incid)))
  # vec.hiv.nnt <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (median(incid.base) - unname(colSums(sim$epi$incid)))
  # vec.gc.nnt <- gc.asympt.tests / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
  # vec.ct.nnt <- ct.asympt.tests / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))
  # vec.syph.nnt <- (syph.asympt.tests) / (median(incid.base.syph) - unname(colSums(sim$epi$incid.syph)))
  # vec.sti.nnt <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (median(incid.base.sti) - unname(colSums(sim$epi$incid.sti)))
  #
  #
  # df$hiv.nnt.low[i] <- quantile(vec.hiv.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$hiv.nnt[i] <- quantile(vec.hiv.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$hiv.nnt.high[i] <- quantile(vec.hiv.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$gc.nnt.low[i] <- quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$gc.nnt[i] <- quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$gc.nnt.high[i] <- quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$ct.nnt.low[i] <- quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$ct.nnt[i] <- quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$ct.nnt.high[i] <- quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$syph.nnt.low[i] <- quantile(vec.syph.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$syph.nnt[i] <- quantile(vec.syph.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$syph.nnt.high[i] <- quantile(vec.syph.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
  #
  # df$sti.nnt.low[i] <- quantile(vec.sti.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
  # df$sti.nnt[i] <- quantile(vec.sti.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
  # df$sti.nnt.high[i] <- quantile(vec.sti.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
#
#   num <- unname(colSums(sim$epi$num.asympt.tx))
#   denom <- unname(colSums(sim$epi$num.asympt.cases))
#   vec.asympt.tx <- num / denom
#   df$asympt.tx.low[i] <- quantile(vec.asympt.tx, probs = qnt.low, na.rm = TRUE, names = FALSE)
#   df$asympt.tx[i] <- quantile(vec.asympt.tx, probs = 0.5, na.rm = TRUE, names = FALSE)
#   df$asympt.tx.high[i] <- quantile(vec.asympt.tx, probs = qnt.high, na.rm = TRUE, names = FALSE)
#
#   num <- unname(colSums(sim$epi$num.rect.tx))
#   denom <- unname(colSums(sim$epi$num.rect.cases))
#   vec.rect.tx <- num / denom
#   df$rect.tx.low[i] <- quantile(vec.rect.tx, probs = qnt.low, na.rm = TRUE, names = FALSE)
#   df$rect.tx[i] <- quantile(vec.rect.tx, probs = 0.5, na.rm = TRUE, names = FALSE)
#   df$rect.tx.high[i] <- quantile(vec.rect.tx, probs = qnt.high, na.rm = TRUE, names = FALSE)

  vec.txasympt <- unname(colSums(sim$epi$txSTI_asympt))
  df$txasympt[i] <- quantile(vec.txasympt, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txasympt.low[i] <- quantile(vec.txasympt, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txasympt.high[i] <- quantile(vec.txasympt, probs = qnt.high, na.rm = TRUE, names = FALSE)

  vec.txCT <- unname(colSums(sim$epi$txCT))
  df$txCT.low[i] <- quantile(vec.txCT, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txCT[i] <- quantile(vec.txCT, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txCT.high[i] <- quantile(vec.txCT, probs = qnt.high, na.rm = TRUE, names = FALSE)

  vec.txGC <- unname(colSums(sim$epi$txGC))
  df$txGC.low[i] <- quantile(vec.txGC, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txGC[i] <- quantile(vec.txGC, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txGC.high[i] <- quantile(vec.txGC, probs = qnt.high, na.rm = TRUE, names = FALSE)

  vec.txsyph <- unname(colSums(sim$epi$txsyph))
  df$txsyph.low[i] <- quantile(vec.txsyph, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txsyph[i] <- quantile(vec.txsyph, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txsyph.high[i] <- quantile(vec.txsyph, probs = qnt.high, na.rm = TRUE, names = FALSE)

  vec.txearlysyph <- unname(colSums(sim$epi$txearlysyph))
  df$txearlysyph.low[i] <- quantile(vec.txearlysyph, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txearlysyph[i] <- quantile(vec.txearlysyph, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txearlysyph.high[i] <- quantile(vec.txearlysyph, probs = qnt.high, na.rm = TRUE, names = FALSE)

  vec.txlatesyph <- unname(colSums(sim$epi$txlatesyph))
  df$txlatesyph.low[i] <- quantile(vec.txlatesyph, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txlatesyph[i] <- quantile(vec.txlatesyph, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txlatesyph.high[i] <- quantile(vec.txlatesyph, probs = qnt.high, na.rm = TRUE, names = FALSE)

  #txCT or txCT_asympt?
  vec.tx.ctpy <- unname(colSums(sim$epi$txCT, na.rm = TRUE) / (sim$epi$num * sim$epi$prev.ct))
  df$txctpy.low[i] <- quantile(vec.tx.ctpy, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txctpy[i] <- quantile(vec.tx.ctpy, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txctpy.high[i] <- quantile(vec.tx.ctpy, probs = qnt.high, na.rm = TRUE, names = FALSE)

  #txGC or txGC_asympt?
  vec.tx.gcpy <- unname(colSums(sim$epi$txGC, na.rm = TRUE) / (sim$epi$num * sim$epi$prev.gc))
  df$txgcpy.low[i] <- quantile(vec.tx.gcpy, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txgcpy[i] <- quantile(vec.tx.gcpy, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txgcpy.high[i] <- quantile(vec.tx.gcpy, probs = qnt.high, na.rm = TRUE, names = FALSE)

  #txsyph or txsyph_asympt or early/late syph?
  vec.tx.syphpy <- unname(colSums(sim$epi$txsyph, na.rm = TRUE) / (sim$epi$num * sim$epi$prev.syph))
  df$txsyphpy.low[i] <- quantile(vec.tx.syphpy, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txsyphpy[i] <- quantile(vec.tx.syphpy, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txsyphpy.high[i] <- quantile(vec.tx.syphpy, probs = qnt.high, na.rm = TRUE, names = FALSE)

  #
  vec.tx.stipy <- unname(colSums(sim$epi$txSTI, na.rm = TRUE) / (sim$epi$num * sim$epi$prev.sti))
  df$txstipy.low[i] <- quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$txstipy[i] <- quantile(vec.tx.stipy, probs = 0.5, na.rm = TRUE, names = FALSE)
  df$txstipy.high[i] <- quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE)

  cat("*")

}

df

write.csv(df, "analysis/STD Table 2.csv")
