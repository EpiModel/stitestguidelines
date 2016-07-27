
# system("scp scripts/burnin/*.burn.[Rs]* hyak:/gscratch/csde/sjenness/sti")
# system("scp scripts/burnin/abc.parms.v1.rda hyak:/gscratch/csde/sjenness/sti")
# system("scp source/*.* hyak:/gscratch/csde/sjenness/sti/source/")

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
sourceDir("source/", FALSE)

## Environmental Arguments
args <- commandArgs(TRUE)
simno <- args[1]
jobno <- args[2]

## Parameters
fsimno <- paste(simno, jobno, sep = ".")

load("est/nwstats.rda")

load("abc.parms.v1.rda")
for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

param <- param_msm(nwstats = st,
                   ai.scale = 1,

                   prep.coverage = 0,
                   prep.start = 2601,

                   rcomp.prob = 0,
                   rcomp.adh.groups = 0:3,
                   rcomp.main.only = FALSE,
                   rcomp.discl.only = FALSE,

                   rgc.tprob = rgc.tprob,
                   ugc.tprob = ugc.tprob,
                   rct.tprob = rct.tprob,
                   uct.tprob = uct.tprob,

                   rgc.sympt.prob = rgc.sympt.prob,
                   ugc.sympt.prob = ugc.sympt.prob,
                   rct.sympt.prob = rct.sympt.prob,
                   uct.sympt.prob = uct.sympt.prob,

                   rgc.dur.asympt = rgc.dur.asympt,
                   ugc.dur.asympt = ugc.dur.asympt,
                   gc.dur.tx = 2,
                   gc.dur.ntx = NULL,

                   rct.dur.asympt = rct.dur.asympt,
                   uct.dur.asympt = uct.dur.asympt,
                   ct.dur.tx = 2,
                   ct.dur.ntx = NULL,

                   gc.prob.cease = 0,
                   ct.prob.cease = 0,

                   gc.sympt.prob.tx = 0.90,
                   ct.sympt.prob.tx = 0.85,
                   gc.asympt.prob.tx = 0,
                   ct.asympt.prob.tx = 0,

                   prep.sti.screen.int = 182,
                   prep.sti.prob.tx = 1,

                   sti.cond.rr = 0.3,

                   hiv.rgc.rr = hiv.rect.rr,
                   hiv.ugc.rr = hiv.ureth.rr,
                   hiv.rct.rr = hiv.rect.rr,
                   hiv.uct.rr = hiv.ureth.rr,
                   hiv.dual.rr = 0)

init <- init_msm(nwstats = st,
                 prev.B = 0.253,
                 prev.W = 0.253,
                 prev.ugc = 0.05,
                 prev.rgc = 0.05,
                 prev.uct = 0.05,
                 prev.rct = 0.05)

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 16,
                       ncores = 16,
                       acts.FUN = acts_sti,
                       condoms.FUN = condoms_sti,
                       initialize.FUN = initialize_sti,
                       prep.FUN = prep_sti,
                       prev.FUN = prevalence_sti,
                       riskhist.FUN = riskhist_sti,
                       position.FUN = position_sti,
                       trans.FUN = trans_sti,
                       stitrans.FUN = sti_trans,
                       stirecov.FUN = sti_recov,
                       stitx.FUN = sti_tx,
                       verbose.FUN = verbose_sti,
                       verbose.int = 100,
                       module.order = c("aging.FUN", "deaths.FUN", "births.FUN",
                                        "test.FUN", "tx.FUN", "prep.FUN",
                                        "progress.FUN", "vl.FUN",
                                        "resim_nets.FUN", "disclose.FUN",
                                        "acts.FUN", "condoms.FUN", "riskhist.FUN",
                                        "position.FUN", "trans.FUN", "stitrans.FUN",
                                        "stirecov.FUN", "stitx.FUN", "prev.FUN"))

## Simulation
netsim_hpc("est/fit.rda", param, init, control, cp.save.int = 1e8,
            save.min = TRUE, save.max = TRUE)

process_simfiles(min.n = 7, outdir = "data/")
