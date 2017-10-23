rm(list=ls())
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
# devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

load("est/nwstats.rda")
load("est/stimod.burnin.rda")
param <- param_msm(nwstats = st,
                   ai.scale = 1.04,

                   # STI acquisition
                   rgc.tprob = 0.4773,
                   ugc.tprob = 0.3819,
                   rct.tprob = 0.2564,
                   uct.tprob = 0.2091,
                   syph.tprob = 0.2533,

                   # HIV acquisition
                   hiv.rgc.rr = 1.75,
                   hiv.ugc.rr = 1.26,
                   hiv.rct.rr = 1.75,
                   hiv.uct.rr = 1.26,
                   hiv.syph.rr = 1.63,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.70,
                   syph.seco.sympt.prob = 0.85,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.80,
                   syph.seco.sympt.prob.tx = 0.80,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   prep.start = 7000,
                   stitest.start = 5201,

                   stitest.elig.model = "sti",

                   prep.coverage = 0.0,
                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.10,
                   stihighrisktest.ct.hivneg.coverage = 0.10,
                   stihighrisktest.syph.hivneg.coverage = 0.10,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.67,
                   stihighrisktest.gc.hivpos.coverage = 0.10,
                   stihighrisktest.ct.hivpos.coverage = 0.10,
                   stihighrisktest.syph.hivpos.coverage = 0.10,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(#simno = 1,
                       start = 5201, nsteps = 5250,
                       nsims = 1, ncores = 1,
                       initialize.FUN = reinit_msm)

sim2 <- netsim(sim, param, init, control)

# netsim_hpc("est/stimod.burnin.rda", param, init, control,
#            compress = TRUE, verbose = FALSE)

process_simfiles(simno = 1, min.n = 4,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 #truncate.at = 5200,
                 vars =
                   c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                     "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                     "ir100.sti",
                     "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                     "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                     "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                     "ir100.syph.tttraj1", "ir100.syph.tttraj2"))


#colMeans(sim2$epi$prop.CT.asympt.tx, na.rm = TRUE)
#colMeans(sim2$epi$prop.GC.asympt.tx, na.rm = TRUE)
#colMeans(sim2$epi$prop.rGC.tx, na.rm = TRUE)
#colMeans(sim2$epi$prop.rCT.tx, na.rm = TRUE)


# Testing/Timing ------------------------------------------------------

#control$bi.mods
#undebug(sti_test_msm)
#undebug(sti_tx_msm)
#undebug(reinit_msm)
#debug(sti_recov)
#undebug(prevalence_msm)

# load("est/stimod.burnin.rda")

at <- 5200
dat <- reinit_msm(sim, param, init, control, s = 1)

at <- at + 1

for (at in 5201:5250) {
  dat <- aging_msm(dat, at)
  dat <- deaths_msm(dat, at)
  dat <- births_msm(dat, at)
  dat <- hiv_test_msm(dat, at)
  dat <- sti_test_msm(dat, at)
  dat <- hiv_tx_msm(dat, at)
  dat <- prep_msm(dat, at)
  dat <- hiv_progress_msm(dat, at)
  dat <- syph_progress_msm(dat, at)
  dat <- hiv_vl_msm(dat, at)
  dat <- simnet_msm(dat, at)
  dat <- hiv_disclose_msm(dat, at)
  dat <- part_msm(dat, at)
  dat <- acts_msm(dat, at)
  dat <- condoms_msm(dat, at)
  dat <- riskhist_prep_msm(dat, at)
  dat <- riskhist_stitest_msm(dat, at)
  dat <- position_msm(dat, at)
  dat <- hiv_trans_msm(dat, at)
  dat <- sti_trans_msm(dat, at)
  dat <- sti_recov_msm(dat, at)
  dat <- sti_tx_msm(dat, at)
  dat <- sti_ept_msm(dat, at)
  dat <- prevalence_msm(dat, at)
  cat("\t", at)
}




#undebug(prep_msm)
#debug(sti_tx)
