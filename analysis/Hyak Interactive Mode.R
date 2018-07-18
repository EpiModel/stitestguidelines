## Burn-in ---------------------------

## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Parameters
load("est/nwstats.rda")
param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   tst.rect.sti.rr = 1,

                   # Correlation
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.52, #0.513, # 0.4773,
                   ugc.tprob = 0.44, # 0.432 # 0.3819,
                   rct.tprob = 0.2797, #0.2794, # 0.2564,
                   uct.tprob = 0.2165, # 0.2161, # 0.2091,
                   syph.tprob = 0.1206, # 0.256

                   # HIV acquisition
                   hiv.rgc.rr = 1.97, #1.75,
                   hiv.ugc.rr = 1.48, #1.27,
                   hiv.rct.rr = 1.97, #1.75,
                   hiv.uct.rr = 1.48, #1.27,
                   hiv.syph.rr = 1.64,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.82,
                   syph.seco.sympt.prob = 0.90,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.85,
                   syph.seco.sympt.prob.tx = 0.85,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0, #0.44, # 0.45
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0, #0.65, #0.67
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 5201,

                   #partlist.start = 1,
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.ugc = 0.002,
                 prev.rgc = 0.002,
                 prev.uct = 0.002,
                 prev.rct = 0.002, # 0.03
                 prev.syph.B = 0, #0.01, # 0.03
                 prev.syph.W = 0) #0.01) # 0.03

control <- control_msm(simno = 1,
                       nsteps = 200,
                       nsims = 3, ncores = 1,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
           save.min = FALSE, save.max = TRUE)

process_simfiles(simno = 1, min.n = 3, compress = TRUE, outdir = "data/")



## Follow-up--------------------------------------------

library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
load("est/nwstats.rda")
# anncov <- 0.1
# hrcov <- 0.1
# anncov <- 0.0
# hrcov <- 0.0
# annint <- 364
# hrint <- 182
# partnercutoff <- 1
# stiasymptx <- 1

# Guidelines
anngcnegcov <- 0.44
anngcposcov <- 0.61
annctnegcov <- 0.44
annctposcov <- 0.61
annsyphnegcov <- 0
annsyphposcov <- 0
hrgcnegcov <- 0
hrgcposcov <- 0
hrctnegcov <- 0
hrctposcov <- 0
hrsyphnegcov <- 0
hrsyphposcov <- 0
annint <- 364
hrint <- 182
partnercutoff <- 1
stiasymptx <- 1

# PAF
rectacq <- 1.97
urethacq <- 1.48
gctrans <- 1.3
cttrans <- 1.3
gccttrans <- 1
allstitrans <- 1

#EPT
eptcov <- 0.1
prov.main.ong <- 0.5
prov.pers.ong <- 0.4
prov.main.end <- 0.4
prov.pers.end <- 0.3
prov.inst <- 0.2
uptake.main <- 0.8
uptake.pers <- 0.8
uptake.inst <- 0.8
eptint <- 60

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   tst.rect.sti.rr = 1,

                   # Correlation
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.5161, #0.513,
                   ugc.tprob = 0.4362, # 0.432
                   rct.tprob = 0.2813, #0.2797,
                   uct.tprob = 0.2195, # 0.2165,
                   syph.tprob = 0, #0.1206,

                   # HIV acquisition
                   hiv.rgc.rr = rectacq,
                   hiv.ugc.rr = urethacq,
                   hiv.rct.rr = rectacq,
                   hiv.uct.rr = urethacq,
                   #hiv.syph.rr = syphacq,

                   # HIV transmission
                   hiv.trans.gc.rr = gctrans,
                   hiv.trans.ct.rr = cttrans,
                   #hiv.trans.syph.rr = syphtrans,
                   hiv.trans.gc.ct.rr = gccttrans,
                   #hiv.trans.gc.syph.rr = gcsyphtrans,
                   #hiv.trans.ct.syph.rr = ctsyphtrans,
                   hiv.trans.allsti.rr = allstitrans,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.82,
                   syph.seco.sympt.prob = 0.90,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.85,
                   syph.seco.sympt.prob.tx = 0.85,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   syph.prim.asympt.prob.tx = stiasymptx,
                   syph.seco.asympt.prob.tx = stiasymptx,
                   syph.earlat.asympt.prob.tx = stiasymptx,
                   syph.latelat.asympt.prob.tx = stiasymptx,
                   syph.tert.asympt.prob.tx = stiasymptx,
                   gc.asympt.prob.tx = stiasymptx,
                   ct.asympt.prob.tx = stiasymptx,

                   partnercut = partnercutoff,

                   stianntest.gc.hivneg.coverage = anngcnegcov,
                   stianntest.ct.hivneg.coverage = annctnegcov,
                   stianntest.syph.hivneg.coverage = annsyphnegcov,
                   stihighrisktest.gc.hivneg.coverage = hrgcnegcov,
                   stihighrisktest.ct.hivneg.coverage = hrctnegcov,
                   stihighrisktest.syph.hivneg.coverage = hrsyphnegcov,
                   stianntest.gc.hivpos.coverage = anngcposcov,
                   stianntest.ct.hivpos.coverage = annctposcov,
                   stianntest.syph.hivpos.coverage = annsyphposcov,
                   stihighrisktest.gc.hivpos.coverage = hrgcposcov,
                   stihighrisktest.ct.hivpos.coverage = hrctposcov,
                   stihighrisktest.syph.hivpos.coverage = hrsyphposcov,
                   prep.coverage = 0,
                   ept.coverage = 0,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 5201,

                   stitest.active.int = annint,
                   sti.highrisktest.int = hrint) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(start = 5201,
                       nsteps = 5720,
                       nsims = 1,
                       ncores = 1,
                       initialize.FUN = reinit_msm,
                       verbose = TRUE)

## Simulation
netsim_hpc("est/stimod.burnin.ept.rda", param, init, control,
           compress = TRUE, verbose = TRUE)

process_simfiles(simno = 60401, min.n = 32,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 truncate.at = 5200,
                 c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                   "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                   "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                   "ir100.sti", "ir100.sti.prep", "ir100.gcct",
                   # "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                   # "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                   # "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                   # "ir100.syph.tttraj1", "ir100.syph.tttraj2",
                   # "ir100.gcct.tttraj1", "ir100.gcct.tttraj2",
                   # "incid.sti.tttraj1", "incid.sti.tttraj2",
                   # "incid.gc.tttraj1", "incid.gc.tttraj2",
                   # "incid.ct.tttraj1", "incid.ct.tttraj2",
                   # "incid.syph.tttraj1", "incid.syph.tttraj2",
                   # "incid.gcct.tttraj1", "incid.gcct.tttraj2",
                   "incid.gc.hivneg", "incid.gc.hivpos",
                   "incid.ct.hivneg", "incid.ct.hivpos",
                   "incid.syph.hivneg", "incid.syph.hivpos",
                   "ir100.gc.hivneg", "ir100.gc.hivpos",
                   "ir100.ct.hivneg", "ir100.ct.hivpos",
                   "ir100.syph.hivneg", "ir100.syph.hivpos",
                   "prop.edges.negneg", "prop.edges.negpos", "prop.edges.pospos",
                   "num.acts.negneg", "num.acts.negpos", "num.acts.pospos",
                   "prop.uai.negneg", "prop.uai.negpos", "prop.uai.pospos",
                   "prop.acts.negneg", "prop.acts.negpos", "prop.acts.pospos",
                   "prop.main.edges.negneg", "prop.main.edges.negpos",
                   "prop.main.edges.pospos", "prop.cas.edges.negneg",
                   "prop.cas.edges.negpos", "prop.cas.edges.pospos",
                   "prop.inst.edges.negneg", "prop.inst.edges.negpos", "prop.inst.edges.pospos",
                   "GCasympttests", "uGCasympttests", "rGCasympttests",
                   "CTasympttests", "uCTasympttests", "rCTasympttests",
                   "syphasympttests", "stiasympttests",
                   "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos",
                   "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos",
                   "syphasympttests.pos", "syphearlyasympttests.pos",
                   "syphlateasympttests.pos", "stiasympttests.pos",
                   # "tt.traj.syph1", "tt.traj.gc1", "tt.traj.ct1",
                   # "tt.traj.syph2", "tt.traj.gc2", "tt.traj.ct2",
                   # "tt.traj.sti1", "tt.traj.sti2",
                   "GCsympttests", "uGCsympttests", "rGCsympttests",
                   "CTsympttests", "uCTsympttests", "rCTsympttests",
                   "syphsympttests", "stisympttests",
                   # "rCTasympttests.tttraj1", "rCTasympttests.tttraj2",
                   # "uCTasympttests.tttraj1", "uCTasympttests.tttraj2",
                   # "CTasympttests.tttraj1", "CTasympttests.tttraj2",
                   # "rGCasympttests.tttraj1", "rGCasympttests.tttraj2",
                   # "uGCasympttests.tttraj1", "uGCasympttests.tttraj2",
                   # "GCasympttests.tttraj1", "GCasympttests.tttraj2",
                   # "syphasympttests.tttraj1", "syphasympttests.tttraj2",
                   # "stiasympttests.tttraj1", "stiasympttests.tttraj2",
                   # "rCTsympttests.tttraj1", "rCTsympttests.tttraj2",
                   # "uCTsympttests.tttraj1", "uCTsympttests.tttraj2",
                   # "CTsympttests.tttraj1", "CTsympttests.tttraj2",
                   # "rGCsympttests.tttraj1", "rGCsympttests.tttraj2",
                   # "uGCsympttests.tttraj1", "uGCsympttests.tttraj2",
                   # "GCsympttests.tttraj1", "GCsympttests.tttraj2",
                   # "syphsympttests.tttraj1", "syphsympttests.tttraj2",
                   # "stisympttests.tttraj1", "stisympttests.tttraj2",
                   # "txGC.tttraj1", "txGC_asympt.tttraj1",
                   # "txGC.tttraj2", "txGC_asympt.tttraj2",
                   # "txCT.tttraj1", "txCT_asympt.tttraj1",
                   # "txCT.tttraj2", "txCT_asympt.tttraj2",
                   # "txsyph.tttraj1", "txsyph_asympt.tttraj1",
                   # "txsyph.tttraj2", "txsyph_asympt.tttraj2",
                   # "txearlysyph.tttraj1", "txlatesyph.tttraj1",
                   # "txearlysyph.tttraj2", "txlatesyph.tttraj2",
                   # "txSTI.tttraj1", "txSTI.tttraj2",
                   "txSTI_asympt.tttraj1", "txSTI_asympt.tttraj2",
                   "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                   'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                   'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                   'test.syph.12mo', 'test.syph.12mo.hivpos', 'test.syph.12mo.hivneg',
                   # "i.prev", "prev.primsecosyph", "prev.syph", "prev.earlysyph", "prev.latesyph",
                   "prev.gc", "prev.rgc", "prev.ugc",
                   "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                   # "prev.dxhiv.dxipssyph", "prev.dxhiv.atdxipssyph",
                   # "prev.primsecosyph.tttraj1", "prev.primsecosyph.tttraj2",
                   # "prev.syph.tttraj1", "prev.syph.tttraj2",
                   # "prev.gcct.tttraj1", "prev.gcct.tttraj2",
                   # "prev.gc.tttraj1", "prev.gc.tttraj2",
                   # "prev.ct.tttraj1", "prev.ct.tttraj2",
                   # "prev.sti.tttraj1", "prev.sti.tttraj2",
                   "prev.rgc.hivneg", "prev.ugc.hivneg",
                   "prev.rct.hivneg", "prev.uct.hivneg",
                   # "prev.primsecosyph.hivneg", "prev.syph.hivneg",
                   "prev.rgc.hivpos","prev.ugc.hivpos",
                   "prev.rct.hivpos", "prev.uct.hivpos",
                   # "prev.primsecosyph.hivpos","prev.syph.hivpos",
                   "prev.hivposmultsti", "prev.hivnegmultsti",
                   #"txearlysyph", "txlatesyph", "txsyph",
                   #"txGC", "txCT",
                   "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                   "tx.gc.prop", "tx.ct.prop", "tx.gcct.prop", #"tx.syph.prop",
                   "rGC_symptstidxtime", "uGC_symptstidxtime", "rCT_symptstidxtime",
                   "uCT_symptstidxtime", #"syph_symptstidxtime",
                   "gc.infect.dur", "ct.infect.dur", "gcct.infect.dur", #"syph.infect.dur",
                   "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                   "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                   "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                   #"cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                   "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                   "deathage", "stiactiveind.prop", "stiactiveind",
                   "recentpartners", "recentpartners.prop",
                   # "num.newearlydiagsyph",
                   # "num.newlatediagsyph", "early.late.syphratio",
                   # "early.late.diagsyphratio",
                   "num.asympt.tx", "num.asympt.cases", "num.rect.cases", "num.rect.tx",
                   "time.hivneg",
                   "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
                   "eptTx", "propindexeptElig", #"eptprop_provided", "eptprop_tx",
                   "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
                   "eptctinfectsti","eptgcinfectundiaghiv", "eptctinfectundiaghiv",
                   "eptgcctinfectundiaghiv",
                   "eptgcinfecthiv", "eptctinfecthiv",
                   "eptgcctinfecthiv",
                   "eptgcctinfecthiv_main", "eptgcctinfecthiv_pers",
                   "eptgcctinfecthiv_inst",
                   "eptgcctinfectundiaghiv_main", "eptgcctinfectundiaghiv_pers",
                   "eptgcctinfectundiaghiv_inst",
                   "eptindexprovided_gc", "eptindexprovided_ct",
                   "eptpartprovided_gc", "eptpartprovided_ct",
                   "eptpartprovided_main", "eptpartprovided_pers",
                   "eptpartprovided_inst", "eptpartuptake_main",
                   "eptpartelig_main", "eptpartelig_pers", "eptpartelig_inst",
                   "eptpartuptake_pers", "eptpartuptake_inst",
                   "eptpartuptake_gc", "eptpartuptake_ct"))#,
                  # "stage.time.ar.ndx","stage.time.ar.dx", "stage.time.af.ndx",
                  # "stage.time.af.dx", "stage.time.early.chronic.ndx",
                  # "stage.time.early.chronic.dx.yrone",
                  # "stage.time.early.chronic.dx.yrstwotolate",
                  # "stage.time.early.chronic.art", "stage.time.late.chronic.ndx",
                  # "stage.time.late.chronic.dx", "stage.time.late.chronic.art",
                  #"stage.time.aids.ndx", "stage.time.aids.dx","stage.time.aids.art"))
