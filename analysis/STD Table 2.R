## STI Testing Guidelines Table 1
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

# Base - No annual or high-risk
#load("data/followup/Guidelines Paper/sim.n3000.rda")
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

# incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
# incid.base.syph.g1 <- unname(colSums(sim.base$epi$incid.syph.tttraj1))
# incid.base.syph.g2 <- unname(colSums(sim.base$epi$incid.syph.tttraj2))
# tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests))
# tests.syph.base.g1 <- unname(colSums(sim.base$epi$syphasympttests.tttraj1))
# tests.syph.base.g2 <- unname(colSums(sim.base$epi$syphasympttests.tttraj2))

incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct))
incid.base.gcct.g1 <- unname(colSums(sim.base$epi$incid.gcct.tttraj1))
incid.base.gcct.g2 <- unname(colSums(sim.base$epi$incid.gcct.tttraj2))
tests.sti.base <- unname(colSums(sim.base$epi$stiasympttests))
tests.sti.base.g1 <- unname(colSums(sim.base$epi$stiasympttests.tttraj1))
tests.sti.base.g2 <- unname(colSums(sim.base$epi$stiasympttests.tttraj2))

# Varying Lower-Risk Coverage
# 3002, 3004, 3006, 3008
# : Annual = 10%, 20%, 30%, 40% increase, 364 days, HR = 0%, 182 days
# Varying Higher-Risk Coverage
# 3018, 3036, 3054, 3072
#  Higher-risk = 0.1 - 1.0 by 0.1, 182 days, Ann = 10%, 364 days

# Varying Lower-Risk Coverage
# 9002, 9004, 9006, 9008
# : Annual = 10%, 20%, 30%, 40% increase, 364 days, HR = 5%, 182 days
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

# syph.incid <- rep(NA, length(sims))
# syph.pia <- rep(NA, length(sims))
# syph.tx.early <- rep(NA, length(sims))
# syph.tx.earlyyr1 <- rep(NA, length(sims))
# syph.tx.late <- rep(NA, length(sims))
# syph.tx.lateyr1 <- rep(NA, length(sims))
# syph.tx <- rep(NA, length(sims))
# syph.txyr1 <- rep(NA, length(sims))
# syph.tx.asympt <- rep(NA, length(sims))
# syph.tx.asymptyr1 <- rep(NA, length(sims))
# syph.tx.sympt <- rep(NA, length(sims))
# syph.tx.symptyr1 <- rep(NA, length(sims))
# syph.txperpy <- rep(NA, length(sims))
# syph.earlytxperpy <- rep(NA, length(sims))
# syph.latetxperpy <- rep(NA, length(sims))
# syph.nnt <- rep(NA, length(sims))
# syph.asympt.tests.py <- rep(NA, length(sims))
# syph.asympt.tests <- rep(NA, length(sims))
#
# syph.incid.g1 <- rep(NA, length(sims))
# syph.pia.g1 <- rep(NA, length(sims))
# syph.asympt.tests.g1 <- rep(NA, length(sims))
# syph.asympt.tests.py.g1 <- rep(NA, length(sims))
# syphtx.g1 <- rep(NA, length(sims))
# syphearlytx.g1 <- rep(NA, length(sims))
# syphlatetx.g1 <- rep(NA, length(sims))
# syphtxpy.g1 <- rep(NA, length(sims))
# syphearlytxpy.g1 <- rep(NA, length(sims))
# syphlatetxpy.g1 <- rep(NA, length(sims))
#
# syph.incid.g2 <- rep(NA, length(sims))
# syph.pia.g2 <- rep(NA, length(sims))
# syph.asympt.tests.g2 <- rep(NA, length(sims))
# syph.asympt.tests.py.g2 <- rep(NA, length(sims))
# syphtx.g2 <- rep(NA, length(sims))
# syphearlytx.g2 <- rep(NA, length(sims))
# syphlatetx.g2 <- rep(NA, length(sims))
# syphtxpy.g2 <- rep(NA, length(sims))
# syphearlytxpy.g2 <- rep(NA, length(sims))
# syphlatetxpy.g2 <- rep(NA, length(sims))
#
# syph.nnt.g1 <- rep(NA, length(sims))
# syph.nnt.g2 <- rep(NA, length(sims))

gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))
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

df <- data.frame(anncov, hrcov, annint, hrint,
                 gc.incid, gc.pia, gc.txyr1, gc.tx, gc.nnt,
                 gc.tx.asymptyr1, gc.tx.asympt, gc.tx.symptyr1, gc.tx.sympt,
                 gctxpy, gc.asympt.tests.py, gc.asympt.tests,
                 ct.incid, ct.pia, ct.txyr1, ct.tx, ct.nnt,
                 ct.tx.asymptyr1, ct.tx.asympt, ct.tx.symptyr1, ct.tx.sympt,
                 cttxpy, ct.asympt.tests.py, ct.asympt.tests,
                 # syph.incid, syph.pia,
                 # syph.tx.earlyyr1, syph.tx.early,
                 # syph.tx.lateyr1, syph.tx.late,
                 # syph.txyr1, syph.tx, syph.nnt,
                 # syph.tx.asymptyr1, syph.tx.asympt,
                 # syph.tx.symptyr1, syph.tx.sympt, syph.txperpy,
                 # syph.earlytxperpy, syph.latetxperpy, syph.asympt.tests.py, syph.asympt.tests,
                 gcct.incid, gcct.pia, sti.txyr1, sti.tx, gcct.nnt,
                 sti.tx.asymptyr1, sti.tx.asympt, sti.tx.symptyr1, sti.tx.sympt,
                 sti.txperpy, sti.asympt.tests.py, sti.asympt.tests,

                 # Group 1
                 gc.incid.g1, gc.pia.g1, gc.nnt.g1,
                 gc.asympt.tests.py.g1, gc.asympt.tests.g1, gctx.g1, gctxpy.g1,
                 ct.incid.g1, ct.pia.g1, ct.nnt.g1,
                 ct.asympt.tests.py.g1, ct.asympt.tests.g1, cttx.g1, cttxpy.g1,
                 # syph.incid.g1, syph.pia.g1, syph.nnt.g1,
                 # syph.asympt.tests.py.g1, syph.asympt.tests.g1, syphtx.g1, syphearlytx.g1,
                 # syphlatetx.g1, syphtxpy.g1, syphearlytxpy.g1, syphlatetxpy.g1,
                 gcct.incid.g1, gcct.pia.g1, gcct.nnt.g1,
                 sti.asympt.tests.py.g1, sti.asympt.tests.g1, tx.sti.g1, txperpy.sti.g1,

                 # Group 2
                 gc.incid.g2, gc.pia.g2, gc.nnt.g2,
                 gc.asympt.tests.py.g2, gc.asympt.tests.g2, gctx.g2, gctxpy.g2,
                 ct.incid.g2, ct.pia.g2, ct.nnt.g2,
                 ct.asympt.tests.py.g2, ct.asympt.tests.g2, cttx.g2, cttxpy.g2,
                 # syph.incid.g2, syph.pia.g2, syph.nnt.g2,
                 # syph.asympt.tests.py.g2, syph.asympt.tests.g2, syphtx.g2,
                 # syphearlytx.g2, syphlatetx.g2, syphtxpy.g2, syphearlytxpy.g2, syphlatetxpy.g2,
                 gcct.incid.g2, gcct.pia.g2, gcct.nnt.g2,
                 sti.asympt.tests.py.g2, sti.asympt.tests.g2, tx.sti.g2, txperpy.sti.g2)

for (i in seq_along(sims)) {

   fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
   load(fn)

   df$anncov[i] <- sim$param$stianntest.ct.hivneg.coverage
   df$hrcov[i] <- sim$param$stihighrisktest.ct.hivpos.coverage
   df$annint[i] <- sim$param$stitest.active.int
   df$hrint[i] <- sim$param$sti.highrisktest.int

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

   # vec.ir.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
   # vec.ir.syph.g1 <- unname(colMeans(tail(sim$epi$ir100.syph.tttraj1, 52)))
   # vec.ir.syph.g2 <- unname(colMeans(tail(sim$epi$ir100.syph.tttraj2, 52)))
   # df$syph.incid[i] <- paste0(round(quantile(vec.ir.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                            " (", round(quantile(vec.ir.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                            " - ", round(quantile(vec.ir.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                            ")")
   # df$syph.incid.g1[i] <- paste0(round(quantile(vec.ir.syph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                               " (", round(quantile(vec.ir.syph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                               " - ", round(quantile(vec.ir.syph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                               ")")
   # df$syph.incid.g2[i] <- paste0(round(quantile(vec.ir.syph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                               " (", round(quantile(vec.ir.syph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                               " - ", round(quantile(vec.ir.syph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                              ")")

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

   #PIA (Cumulative)
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

   # incid.syph <- unname(colSums(sim$epi$incid.syph))
   # vec.nia.syph <- incid.base.syph - incid.syph
   # vec.pia.syph <- vec.nia.syph/incid.base.syph
   #
   # incid.syph.g1 <- unname(colSums(sim$epi$incid.syph.tttraj1))
   # vec.nia.syph.g1 <- round(incid.base.syph.g1 - incid.syph.g1, 1)
   # vec.pia.syph.g1 <- vec.nia.syph.g1/incid.syph.g1
   #
   # incid.syph.g2 <- unname(colSums(sim$epi$incid.syph.tttraj2))
   # vec.nia.syph.g2 <- round(incid.base.syph.g2 - incid.syph.g2, 1)
   # vec.pia.syph.g2 <- vec.nia.syph.g2/incid.syph.g2

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
   # df$syph.pia[i] <- paste0(round(quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                          " (", round(quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                          " - ", round(quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                          ")")
   #
   # df$syph.pia.g1[i] <- paste0(round(quantile(vec.pia.syph.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                             " (", round(quantile(vec.pia.syph.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                             " - ", round(quantile(vec.pia.syph.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                             ")")
   # df$syph.pia.g2[i] <- paste0(round(quantile(vec.pia.syph.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                             " (", round(quantile(vec.pia.syph.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                             " - ", round(quantile(vec.pia.syph.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                             ")")

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

   # vec.txearlysyph <- unname(colSums(sim$epi$txearlysyph))
   # vec.txearlysyphyr1 <- as.numeric(colSums(head(sim$epi$txearlysyph, 52)))
   # df$syph.tx.early[i] <- paste0(round(quantile(vec.txearlysyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                             " (", round(quantile(vec.txearlysyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                             " - ", round(quantile(vec.txearlysyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                             ")")
   # df$syph.tx.earlyyr1[i] <- paste0(round(quantile(vec.txearlysyphyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                               " (", round(quantile(vec.txearlysyphyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                               " - ", round(quantile(vec.txearlysyphyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                               ")")
   #
   # vec.txlatesyph <- unname(colSums(sim$epi$txlatesyph))
   # vec.txlatesyphyr1 <- as.numeric(colSums(head(sim$epi$txlatesyph, 52)))
   # df$syph.tx.late[i] <- paste0(round(quantile(vec.txlatesyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                            " (", round(quantile(vec.txlatesyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                            " - ", round(quantile(vec.txlatesyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                            ")")
   # df$syph.tx.lateyr1[i] <- paste0(round(quantile(vec.txlatesyphyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                              " (", round(quantile(vec.txlatesyphyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                              " - ", round(quantile(vec.txlatesyphyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                              ")")
   #
   # vec.txsyph <- unname(colSums(sim$epi$txsyph))
   # vec.txsyphyr1 <- as.numeric(colSums(head(sim$epi$txearlysyph, 52)))
   # df$syph.txyr1[i] <- paste0(round(quantile(vec.txsyphyr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                         " (", round(quantile(vec.txsyphyr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                         " - ", round(quantile(vec.txsyphyr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                         ")")
   # df$syph.tx[i] <- paste0(round(quantile(vec.txsyph, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                               " (", round(quantile(vec.txsyph, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                               " - ", round(quantile(vec.txsyph, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                               ")")
   #
   # vec.txsyph.asympt <- unname(colSums(sim$epi$txsyph_asympt))
   # vec.txsyph.asympt.yr1 <- as.numeric(colSums(head(sim$epi$txsyph_asympt, 52)))
   # vec.txsyph.sympt <- unname(colSums(sim$epi$txsyph - sim$epi$txsyph_asympt))
   # vec.txsyph.sympt.yr1 <- as.numeric(colSums(head(sim$epi$txsyph - sim$epi$txsyph_asympt, 52)))
   #
   # df$syph.tx.asympt[i] <- paste0(round(quantile(vec.txsyph.asympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                              " (", round(quantile(vec.txsyph.asympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                              " - ", round(quantile(vec.txsyph.asympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                              ")")
   # df$syph.tx.asymptyr1[i] <- paste0(round(quantile(vec.txsyph.asympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                                 " (", round(quantile(vec.txsyph.asympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                                 " - ", round(quantile(vec.txsyph.asympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                                 ")")
   # df$syph.tx.sympt[i] <- paste0(round(quantile(vec.txsyph.sympt, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                             " (", round(quantile(vec.txsyph.sympt, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                             " - ", round(quantile(vec.txsyph.sympt, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                             ")")
   # df$syph.tx.symptyr1[i] <- paste0(round(quantile(vec.txsyph.sympt.yr1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                                " (", round(quantile(vec.txsyph.sympt.yr1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                                " - ", round(quantile(vec.txsyph.sympt.yr1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                                ")")


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

   # syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
   # syph.sympt.tests <- unname(colSums(sim$epi$syphsympttests, na.rm = TRUE))
   # syph.tests <- syph.asympt.tests + syph.sympt.tests
   #
   # syph.asympt.tests.g1 <- unname(colSums(sim$epi$syphasympttests.tttraj1, na.rm = TRUE))
   # syph.sympt.tests.g1 <- unname(colSums(sim$epi$syphsympttests.tttraj1, na.rm = TRUE))
   # syph.tests.g1 <- syph.asympt.tests.g1 + syph.sympt.tests.g1
   #
   # syph.asympt.tests.g2 <- unname(colSums(sim$epi$syphasympttests.tttraj2, na.rm = TRUE))
   # syph.sympt.tests.g2 <- unname(colSums(sim$epi$syphsympttests.tttraj2, na.rm = TRUE))
   # syph.tests.g2 <- syph.asympt.tests.g2 + syph.sympt.tests.g2

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
   # syph.asympt.tests.py <-  52 * (syph.asympt.tests / py)
   # syph.asympt.tests.py.g1 <-  52 * (syph.asympt.tests.g1 / py.g1)
   # syph.asympt.tests.py.g2 <-  52 * (syph.asympt.tests.g2 / py.g2)
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

   # df$syph.asympt.tests[i] <- paste0(round(quantile(syph.asympt.tests, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                                   " (", round(quantile(syph.asympt.tests, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                                   " - ", round(quantile(syph.asympt.tests, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                                   ")")
   # df$syph.asympt.tests.g1[i] <- paste0(round(quantile(syph.asympt.tests.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                                      " (", round(quantile(syph.asympt.tests.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                                      " - ", round(quantile(syph.asympt.tests.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                                      ")")
   # df$syph.asympt.tests.g2[i] <- paste0(round(quantile(syph.asympt.tests.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 0),
   #                                      " (", round(quantile(syph.asympt.tests.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 0),
   #                                      " - ", round(quantile(syph.asympt.tests.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 0),
   #                                      ")")
   # df$syph.asympt.tests.py[i] <- paste0(round(quantile(syph.asympt.tests.py, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                                      " (", round(quantile(syph.asympt.tests.py, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                                      " - ", round(quantile(syph.asympt.tests.py, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                                      ")")
   #
   # df$syph.asympt.tests.py.g1[i] <- paste0(round(quantile(syph.asympt.tests.py.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                                         " (", round(quantile(syph.asympt.tests.py.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                                         " - ", round(quantile(syph.asympt.tests.py.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                                         ")")
   # df$syph.asympt.tests.py.g2[i] <- paste0(round(quantile(syph.asympt.tests.py.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                                         " (", round(quantile(syph.asympt.tests.py.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                                         " - ", round(quantile(syph.asympt.tests.py.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                                         ")")


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
   df$ct.txperpy[i] <- paste0(round(quantile(vec.tx.ctpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
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
   # vec.tx.syphpy <- unname(colMeans(52 * sim$epi$txsyph / (sim$epi$num * sim$epi$prev.syph)))
   # df$syph.txperpy[i] <- paste0(round(quantile(vec.tx.syphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                             " (", round(quantile(vec.tx.syphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                             " - ", round(quantile(vec.tx.syphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                             ")")
   #
   #
   #  # Prev.primsecosyph (Incub, prim, syph)
   # df.prev.earlysyph <- sim$epi$txearlysyph[1:521,] / (sim$epi$num[1:521,] * sim$epi$prev.primsecosyph[1:521,])
   #
   # # Remove NaNs
   # for (j in 1:ncol(df.prev.earlysyph)) {
   #
   #   df.prev.earlysyph[which(is.nan(df.prev.earlysyph[, j])), j] <- 0.0
   #
   # }
   # vec.tx.earlysyphpy <- unname(colMeans(52 * df.prev.earlysyph))
   #
   # #vec.tx.earlysyphpy <- unname(colMeans(52 * sim$epi$txearlysyph / (sim$epi$num * sim$epi$prev.primsecosyph)))
   # df$syph.earlytxperpy[i] <- paste0(round(quantile(vec.tx.earlysyphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                                  " (", round(quantile(vec.tx.earlysyphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                                  " - ", round(quantile(vec.tx.earlysyphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                                  ")")
   #
   # # Late syphilis prevalence = Prev.syph - Prev.primsecosyph
   # df.prev.latesyph <- sim$epi$txlatesyph[1:521,] / (sim$epi$num[1:521,] * (sim$epi$prev.syph[1:521, ] - sim$epi$prev.primsecosyph[1:521,]))
   #
   # # Remove NaNs
   # for (j in 1:ncol(df.prev.latesyph)) {
   #
   #   df.prev.latesyph[which(is.nan(df.prev.latesyph[, j])), j] <- 0.0
   #
   # }
   # vec.tx.latesyphpy <- unname(colMeans(52 * df.prev.latesyph))
   # #vec.tx.latesyphpy <- unname(colMeans(52 * sim$epi$txlatesyph / (sim$epi$num * (sim$epi$prev.syph - sim$epi$prev.primsecosyph))))
   # df$syph.latetxperpy[i] <- paste0(round(quantile(vec.tx.latesyphpy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                                  " (", round(quantile(vec.tx.latesyphpy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                                  " - ", round(quantile(vec.tx.latesyphpy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                                  ")")

   vec.tx.stipy <- unname(colMeans(52 * sim$epi$txSTI / (sim$epi$num * sim$epi$prev.sti)))
   df$sti.txperpy[i] <- paste0(round(quantile(vec.tx.stipy, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.tx.stipy, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.tx.stipy, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

   # Number needed to treat
   # incremental tests / number of infections averted
   vec.gc.nnt <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - incid.gc)
   vec.gc.nnt.g1 <- (gc.asympt.tests.g1 - tests.gc.base.g1) / (incid.base.gc.g1 - incid.gc.g1)
   vec.gc.nnt.g2 <- (gc.asympt.tests.g2 - tests.gc.base.g2) / (incid.base.gc.g2 - incid.gc.g2)

   vec.ct.nnt <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - incid.ct)
   vec.ct.nnt.g1 <- (ct.asympt.tests.g1 - tests.ct.base.g1) / (incid.base.ct.g1 - incid.ct.g1)
   vec.ct.nnt.g2 <- (ct.asympt.tests.g2 - tests.ct.base.g2) / (incid.base.ct.g2 - incid.ct.g2)

   # vec.syph.nnt <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - incid.syph)
   # vec.syph.nnt.g1 <- (syph.asympt.tests.g1  - tests.syph.base.g1) / (incid.base.syph.g1 - incid.syph.g1)
   # vec.syph.nnt.g2 <- (syph.asympt.tests.g2  - tests.syph.base.g2) / (incid.base.syph.g2 - incid.syph.g2)

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

   # df$syph.nnt[i] <- paste0(round(quantile(vec.syph.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                          " (", round(quantile(vec.syph.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                          " - ", round(quantile(vec.syph.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                          ")")
   # df$syph.nnt.g1[i] <- paste0(round(quantile(vec.syph.nnt.g1, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                             " (", round(quantile(vec.syph.nnt.g1, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                             " - ", round(quantile(vec.syph.nnt.g1, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                             ")")
   # df$syph.nnt.g2[i] <- paste0(round(quantile(vec.syph.nnt.g2, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
   #                             " (", round(quantile(vec.syph.nnt.g2, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
   #                             " - ", round(quantile(vec.syph.nnt.g2, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
   #                             ")")

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

write.csv(df, "analysis/STD Table 2.csv")
