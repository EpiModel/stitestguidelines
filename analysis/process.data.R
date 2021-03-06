## Process STI Testing Guidelines Data

# Bulk on Hyak -----------------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
#source("analysis/fx.R")

# ( fn <- list.files("data/followup", full.names = TRUE) )
# # fn <- list.files(pattern = "data/followup/n[3-4][0-9][0-9][0-9].rda")
# for (i in fn) {
#     load(i)
#     sim <- truncate_sim(sim, at = 2600)
#     save(sim, file = i, compress = TRUE)
#     cat("*")
# }

# truncate and limit sim 3000/4000 files on Hyak
# cd /gscratch/csde/kweiss2/sti/data
# module load r_3.2.4
# R
fn <- list.files(pattern = "n[8-9][0-9][0-9][0-9].rda")
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 5200)
  vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                   "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                   "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                   "ir100.sti", "ir100.sti.prep", "ir100.gcct",
                   "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                   "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                   "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                   "ir100.syph.tttraj1", "ir100.syph.tttraj2",
                   "ir100.gcct.tttraj1", "ir100.gcct.tttraj2",
                   "incid.sti.tttraj1", "incid.sti.tttraj2",
                   "incid.gc.tttraj1", "incid.gc.tttraj2",
                   "incid.ct.tttraj1", "incid.ct.tttraj2",
                   "incid.syph.tttraj1", "incid.syph.tttraj2",
                   "incid.gcct.tttraj1", "incid.gcct.tttraj2",
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
                   "tt.traj.syph1", "tt.traj.gc1", "tt.traj.ct1",
                   "tt.traj.syph2", "tt.traj.gc2", "tt.traj.ct2",
                   "tt.traj.sti1", "tt.traj.sti2",
                   "GCsympttests", "uGCsympttests", "rGCsympttests",
                   "CTsympttests", "uCTsympttests", "rCTsympttests",
                   "syphsympttests", "stisympttests",
                   "rCTasympttests.tttraj1", "rCTasympttests.tttraj2",
                   "uCTasympttests.tttraj1", "uCTasympttests.tttraj2",
                   "CTasympttests.tttraj1", "CTasympttests.tttraj2",
                   "rGCasympttests.tttraj1", "rGCasympttests.tttraj2",
                   "uGCasympttests.tttraj1", "uGCasympttests.tttraj2",
                   "GCasympttests.tttraj1", "GCasympttests.tttraj2",
                   "syphasympttests.tttraj1", "syphasympttests.tttraj2",
                   "stiasympttests.tttraj1", "stiasympttests.tttraj2",
                   "rCTsympttests.tttraj1", "rCTsympttests.tttraj2",
                   "uCTsympttests.tttraj1", "uCTsympttests.tttraj2",
                   "CTsympttests.tttraj1", "CTsympttests.tttraj2",
                   "rGCsympttests.tttraj1", "rGCsympttests.tttraj2",
                   "uGCsympttests.tttraj1", "uGCsympttests.tttraj2",
                   "GCsympttests.tttraj1", "GCsympttests.tttraj2",
                   "syphsympttests.tttraj1", "syphsympttests.tttraj2",
                   "stisympttests.tttraj1", "stisympttests.tttraj2",
                   "txGC.tttraj1", "txGC_asympt.tttraj1",
                   "txGC.tttraj2", "txGC_asympt.tttraj2",
                   "txCT.tttraj1", "txCT_asympt.tttraj1",
                   "txCT.tttraj2", "txCT_asympt.tttraj2",
                   "txsyph.tttraj1", "txsyph_asympt.tttraj1",
                   "txsyph.tttraj2", "txsyph_asympt.tttraj2",
                   "txearlysyph.tttraj1", "txlatesyph.tttraj1",
                   "txearlysyph.tttraj2", "txlatesyph.tttraj2",
                   "txSTI.tttraj1", "txSTI.tttraj2",
                   "txSTI_asympt.tttraj1", "txSTI_asympt.tttraj2",
                   "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                   'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                   'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                   'test.syph.12mo', 'test.syph.12mo.hivpos', 'test.syph.12mo.hivneg',
                   "i.prev", "prev.primsecosyph", "prev.syph", "prev.earlysyph", "prev.latesyph",
                   "prev.gc", "prev.rgc", "prev.ugc",
                   "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                   "prev.dxhiv.dxipssyph", "prev.dxhiv.atdxipssyph",
                   "prev.primsecosyph.tttraj1", "prev.primsecosyph.tttraj2",
                   "prev.syph.tttraj1", "prev.syph.tttraj2",
                   "prev.gcct.tttraj1", "prev.gcct.tttraj2",
                   "prev.gc.tttraj1", "prev.gc.tttraj2",
                   "prev.ct.tttraj1", "prev.ct.tttraj2",
                   "prev.sti.tttraj1", "prev.sti.tttraj2",
                   "prev.rgc.hivneg", "prev.ugc.hivneg",
                   "prev.rct.hivneg", "prev.uct.hivneg",
                   "prev.primsecosyph.hivneg", "prev.syph.hivneg",
                   "prev.rgc.hivpos","prev.ugc.hivpos",
                   "prev.rct.hivpos", "prev.uct.hivpos",
                   "prev.primsecosyph.hivpos","prev.syph.hivpos",
                   "prev.hivposmultsti", "prev.hivnegmultsti",
                   "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                   "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                   "tx.gc.prop", "tx.ct.prop", "tx.gcct.prop", "tx.syph.prop",
                   "rGC_symptstidxtime", "uGC_symptstidxtime", "rCT_symptstidxtime",
                   "uCT_symptstidxtime", "syph_symptstidxtime",
                   "gc.infect.dur", "ct.infect.dur", "gcct.infect.dur", "syph.infect.dur",
                   "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                   "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                   "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                   "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                   "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                   "deathage", "stiactiveind.prop", "stiactiveind",
                   "recentpartners", "recentpartners.prop",
                   "num.newearlydiagsyph",
                   "num.newlatediagsyph", "early.late.syphratio",
                   "early.late.diagsyphratio",
                   "num.asympt.tx", "num.asympt.cases", "num.rect.cases", "num.rect.tx",
                   "time.hivneg",
                   "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
                   "eptTx", "propindexeptElig", "eptprop_provided", "eptprop_tx",
                   "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
                   "eptctinfectsti","eptgcinfecthiv", "eptctinfecthiv",
                   "stage.time.ar.ndx","stage.time.ar.dx", "stage.time.af.ndx",
                   "stage.time.af.dx", "stage.time.early.chronic.ndx",
                   "stage.time.early.chronic.dx.yrone",
                   "stage.time.early.chronic.dx.yrstwotolate",
                   "stage.time.early.chronic.art", "stage.time.late.chronic.ndx",
                   "stage.time.late.chronic.dx", "stage.time.late.chronic.art",
                   "stage.time.aids.ndx", "stage.time.aids.dx","stage.time.aids.art")

  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  out.fn <- paste0("followup/", i)
  save(sim, file = out.fn, compress = "gzip")
  file.remove(i)
  cat(i, "\n")
}


### 1 by 1 processing on Hyak - cd /gscratch/csde/kweiss2/sti/data -------------
rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
load("sim.n3037.rda")
sim <- truncate_sim(sim, at = 5200)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.sti.prep", "ir100.gcct",
                 "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                 "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                 "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                 "ir100.syph.tttraj1", "ir100.syph.tttraj2",
                 "ir100.gcct.tttraj1", "ir100.gcct.tttraj2",
                 "incid.sti.tttraj1", "incid.sti.tttraj2",
                 "incid.gc.tttraj1", "incid.gc.tttraj2",
                 "incid.ct.tttraj1", "incid.ct.tttraj2",
                 "incid.syph.tttraj1", "incid.syph.tttraj2",
                 "incid.gcct.tttraj1", "incid.gcct.tttraj2",
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
                 "tt.traj.syph1", "tt.traj.gc1", "tt.traj.ct1",
                 "tt.traj.syph2", "tt.traj.gc2", "tt.traj.ct2",
                 "tt.traj.sti1", "tt.traj.sti2",
                 "GCsympttests", "uGCsympttests", "rGCsympttests",
                 "CTsympttests", "uCTsympttests", "rCTsympttests",
                 "syphsympttests", "stisympttests",
                 "rCTasympttests.tttraj1", "rCTasympttests.tttraj2",
                 "uCTasympttests.tttraj1", "uCTasympttests.tttraj2",
                 "CTasympttests.tttraj1", "CTasympttests.tttraj2",
                 "rGCasympttests.tttraj1", "rGCasympttests.tttraj2",
                 "uGCasympttests.tttraj1", "uGCasympttests.tttraj2",
                 "GCasympttests.tttraj1", "GCasympttests.tttraj2",
                 "syphasympttests.tttraj1", "syphasympttests.tttraj2",
                 "stiasympttests.tttraj1", "stiasympttests.tttraj2",
                 "rCTsympttests.tttraj1", "rCTsympttests.tttraj2",
                 "uCTsympttests.tttraj1", "uCTsympttests.tttraj2",
                 "CTsympttests.tttraj1", "CTsympttests.tttraj2",
                 "rGCsympttests.tttraj1", "rGCsympttests.tttraj2",
                 "uGCsympttests.tttraj1", "uGCsympttests.tttraj2",
                 "GCsympttests.tttraj1", "GCsympttests.tttraj2",
                 "syphsympttests.tttraj1", "syphsympttests.tttraj2",
                 "stisympttests.tttraj1", "stisympttests.tttraj2",
                 "txGC.tttraj1", "txGC_asympt.tttraj1",
                 "txGC.tttraj2", "txGC_asympt.tttraj2",
                 "txCT.tttraj1", "txCT_asympt.tttraj1",
                 "txCT.tttraj2", "txCT_asympt.tttraj2",
                 "txsyph.tttraj1", "txsyph_asympt.tttraj1",
                 "txsyph.tttraj2", "txsyph_asympt.tttraj2",
                 "txearlysyph.tttraj1", "txlatesyph.tttraj1",
                 "txearlysyph.tttraj2", "txlatesyph.tttraj2",
                 "txSTI.tttraj1", "txSTI.tttraj2",
                 "txSTI_asympt.tttraj1", "txSTI_asympt.tttraj2",
                 "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                 'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                 'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                 'test.syph.12mo', 'test.syph.12mo.hivpos', 'test.syph.12mo.hivneg',
                 "i.prev", "prev.primsecosyph", "prev.syph", "prev.earlysyph", "prev.latesyph",
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                 "prev.dxhiv.dxipssyph", "prev.dxhiv.atdxipssyph",
                 "prev.primsecosyph.tttraj1", "prev.primsecosyph.tttraj2",
                 "prev.syph.tttraj1", "prev.syph.tttraj2",
                 "prev.gcct.tttraj1", "prev.gcct.tttraj2",
                 "prev.gc.tttraj1", "prev.gc.tttraj2",
                 "prev.ct.tttraj1", "prev.ct.tttraj2",
                 "prev.sti.tttraj1", "prev.sti.tttraj2",
                 "prev.rgc.hivneg", "prev.ugc.hivneg",
                 "prev.rct.hivneg", "prev.uct.hivneg",
                 "prev.primsecosyph.hivneg", "prev.syph.hivneg",
                 "prev.rgc.hivpos","prev.ugc.hivpos",
                 "prev.rct.hivpos", "prev.uct.hivpos",
                 "prev.primsecosyph.hivpos","prev.syph.hivpos",
                 "prev.hivposmultsti", "prev.hivnegmultsti",
                 "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                 "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                 "tx.gc.prop", "tx.ct.prop", "tx.gcct.prop", "tx.syph.prop",
                 "rGC_symptstidxtime", "uGC_symptstidxtime", "rCT_symptstidxtime",
                 "uCT_symptstidxtime", "syph_symptstidxtime",
                 "gc.infect.dur", "ct.infect.dur", "gcct.infect.dur", "syph.infect.dur",
                 "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                 "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                 "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                 "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                 "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                 "deathage", "stiactiveind.prop", "stiactiveind",
                 "recentpartners", "recentpartners.prop",
                 "num.newearlydiagsyph",
                 "num.newlatediagsyph", "early.late.syphratio",
                 "early.late.diagsyphratio",
                 "num.asympt.tx", "num.asympt.cases", "num.rect.cases", "num.rect.tx",
                 "time.hivneg",
                 "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
                 "eptTx", "propindexeptElig", "eptprop_provided", "eptprop_tx",
                 "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
                 "eptctinfectsti","eptgcinfecthiv", "eptctinfecthiv",
                 "stage.time.ar.ndx","stage.time.ar.dx", "stage.time.af.ndx",
                 "stage.time.af.dx", "stage.time.early.chronic.ndx",
                 "stage.time.early.chronic.dx.yrone",
                 "stage.time.early.chronic.dx.yrstwotolate",
                 "stage.time.early.chronic.art", "stage.time.late.chronic.ndx",
                 "stage.time.late.chronic.dx", "stage.time.late.chronic.art",
                 "stage.time.aids.ndx", "stage.time.aids.dx","stage.time.aids.art")
i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "followup/sim.n3037.rda", compress = "gzip")


## Locally merge files --------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

sim <- merge_simfiles(simno = 3118, indir = "data/", ftype = "min")
sim <- truncate_sim(sim, at = 5200)
vars.needed <-  c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                  "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                  "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                  "ir100.sti", "ir100.sti.prep", "ir100.gcct",
                  "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                  "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                  "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                  "ir100.syph.tttraj1", "ir100.syph.tttraj2",
                  "ir100.gcct.tttraj1", "ir100.gcct.tttraj2",
                  "incid.sti.tttraj1", "incid.sti.tttraj2",
                  "incid.gc.tttraj1", "incid.gc.tttraj2",
                  "incid.ct.tttraj1", "incid.ct.tttraj2",
                  "incid.syph.tttraj1", "incid.syph.tttraj2",
                  "incid.gcct.tttraj1", "incid.gcct.tttraj2",
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
                  "tt.traj.syph1", "tt.traj.gc1", "tt.traj.ct1",
                  "tt.traj.syph2", "tt.traj.gc2", "tt.traj.ct2",
                  "tt.traj.sti1", "tt.traj.sti2",
                  "GCsympttests", "uGCsympttests", "rGCsympttests",
                  "CTsympttests", "uCTsympttests", "rCTsympttests",
                  "syphsympttests", "stisympttests",
                  "rCTasympttests.tttraj1", "rCTasympttests.tttraj2",
                  "uCTasympttests.tttraj1", "uCTasympttests.tttraj2",
                  "CTasympttests.tttraj1", "CTasympttests.tttraj2",
                  "rGCasympttests.tttraj1", "rGCasympttests.tttraj2",
                  "uGCasympttests.tttraj1", "uGCasympttests.tttraj2",
                  "GCasympttests.tttraj1", "GCasympttests.tttraj2",
                  "syphasympttests.tttraj1", "syphasympttests.tttraj2",
                  "stiasympttests.tttraj1", "stiasympttests.tttraj2",
                  "rCTsympttests.tttraj1", "rCTsympttests.tttraj2",
                  "uCTsympttests.tttraj1", "uCTsympttests.tttraj2",
                  "CTsympttests.tttraj1", "CTsympttests.tttraj2",
                  "rGCsympttests.tttraj1", "rGCsympttests.tttraj2",
                  "uGCsympttests.tttraj1", "uGCsympttests.tttraj2",
                  "GCsympttests.tttraj1", "GCsympttests.tttraj2",
                  "syphsympttests.tttraj1", "syphsympttests.tttraj2",
                  "stisympttests.tttraj1", "stisympttests.tttraj2",
                  "txGC.tttraj1", "txGC_asympt.tttraj1",
                  "txGC.tttraj2", "txGC_asympt.tttraj2",
                  "txCT.tttraj1", "txCT_asympt.tttraj1",
                  "txCT.tttraj2", "txCT_asympt.tttraj2",
                  "txsyph.tttraj1", "txsyph_asympt.tttraj1",
                  "txsyph.tttraj2", "txsyph_asympt.tttraj2",
                  "txearlysyph.tttraj1", "txlatesyph.tttraj1",
                  "txearlysyph.tttraj2", "txlatesyph.tttraj2",
                  "txSTI.tttraj1", "txSTI.tttraj2",
                  "txSTI_asympt.tttraj1", "txSTI_asympt.tttraj2",
                  "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                  'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                  'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                  'test.syph.12mo', 'test.syph.12mo.hivpos', 'test.syph.12mo.hivneg',
                  "i.prev", "prev.primsecosyph", "prev.syph", "prev.earlysyph", "prev.latesyph",
                  "prev.gc", "prev.rgc", "prev.ugc",
                  "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                  "prev.dxhiv.dxipssyph", "prev.dxhiv.atdxipssyph",
                  "prev.primsecosyph.tttraj1", "prev.primsecosyph.tttraj2",
                  "prev.syph.tttraj1", "prev.syph.tttraj2",
                  "prev.gcct.tttraj1", "prev.gcct.tttraj2",
                  "prev.gc.tttraj1", "prev.gc.tttraj2",
                  "prev.ct.tttraj1", "prev.ct.tttraj2",
                  "prev.sti.tttraj1", "prev.sti.tttraj2",
                  "prev.rgc.hivneg", "prev.ugc.hivneg",
                  "prev.rct.hivneg", "prev.uct.hivneg",
                  "prev.primsecosyph.hivneg", "prev.syph.hivneg",
                  "prev.rgc.hivpos","prev.ugc.hivpos",
                  "prev.rct.hivpos", "prev.uct.hivpos",
                  "prev.primsecosyph.hivpos","prev.syph.hivpos",
                  "prev.hivposmultsti", "prev.hivnegmultsti",
                  "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                  "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                  "tx.gc.prop", "tx.ct.prop", "tx.gcct.prop", "tx.syph.prop",
                  "rGC_symptstidxtime", "uGC_symptstidxtime", "rCT_symptstidxtime",
                  "uCT_symptstidxtime", "syph_symptstidxtime",
                  "gc.infect.dur", "ct.infect.dur", "gcct.infect.dur", "syph.infect.dur",
                  "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                  "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                  "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                  "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                  "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                  "deathage", "stiactiveind.prop", "stiactiveind",
                  "recentpartners", "recentpartners.prop",
                  "num.newearlydiagsyph",
                  "num.newlatediagsyph", "early.late.syphratio",
                  "early.late.diagsyphratio",
                  "num.asympt.tx", "num.asympt.cases", "num.rect.cases", "num.rect.tx",
                  "time.hivneg",
                  "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
                  "eptTx", "propindexeptElig", "eptprop_provided", "eptprop_tx",
                  "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
                  "eptctinfectsti","eptgcinfecthiv", "eptctinfecthiv",
                  "stage.time.ar.ndx","stage.time.ar.dx", "stage.time.af.ndx",
                  "stage.time.af.dx", "stage.time.early.chronic.ndx",
                  "stage.time.early.chronic.dx.yrone",
                  "stage.time.early.chronic.dx.yrstwotolate",
                  "stage.time.early.chronic.art", "stage.time.late.chronic.ndx",
                  "stage.time.late.chronic.dx", "stage.time.late.chronic.art",
                  "stage.time.aids.ndx", "stage.time.aids.dx","stage.time.aids.art")

i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "data/followup/sim.3118.rda", compress = "gzip")


# 1 by 1 burnin on Hyak---------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
sim <- merge_simfiles(simno = 1000, indir = "data/", ftype = "max")
save(sim, file = "data/sim.n1000.rda", compress = "gzip")


#### Merge 1 by 1 on Hyak ------------------------------------------------------
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

sim <- merge_simfiles(simno = 6041, indir = "data/", ftype = "min")
sim <- truncate_sim(sim, at = 5200)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.gcct",
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
                 "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                 'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                 'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                 "prev.rgc.hivneg", "prev.ugc.hivneg",
                 "prev.rct.hivneg", "prev.uct.hivneg",
                 "prev.rgc.hivpos","prev.ugc.hivpos",
                 "prev.rct.hivpos", "prev.uct.hivpos",
                 "prev.hivposmultsti", "prev.hivnegmultsti",
                 "txGC", "txCT",
                 "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                 "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                 "cell1_rectureth", "cell2_rectureth", "cell3_rectureth", "cell4_rectureth",
                 "cell1_newinf", "cell2_newinf", "cell3_newinf", "cell4_newinf",
                 "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                 "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                 "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                 "stiactiveind.prop", "stiactiveind",
                 "recentpartners", "recentpartners.prop",
                 "rGC_hiv", "uGC_hiv", "rCT_hiv", "uCT_hiv", "rGC_hiv_acq",
                 "uGC_hiv_acq", "rCT_hiv_acq", "uCT_hiv_acq", "rGC_hiv_trans",
                 "uGC_hiv_trans", "rCT_hiv_trans", "uCT_hiv_trans",
                 "rGC_hiv_trans_events_perperson",
                 "uGC_hiv_trans_events_perperson",
                 "rCT_hiv_trans_events_perperson",
                 "uCT_hiv_trans_events_perperson")
i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "data/sim.n3000.rda", compress = "gzip")

### More Hyak merge------------------------------------------------------------------

# sims <- c(3001:3009, 3018, 3027, 3036, 3045, 3054, 3063, 3072, 3081, 3090, 3099, 3108, 3117, 3126, 3135, 3144, 3153, 3162, 3171, 3180, 3189:3198, 3221:3513)
#sims <- c(4001:4009, 4018, 4027, 4036, 4045, 4054, 4063, 4072, 4081, 4090, 4099, 4108, 4117, 4126, 4135, 4144, 4153, 4162, 4171, 4180, 6001:6009)
rm(list = ls())
library("EpiModel")
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

sims <- c(9000:9416)
for (i in sims) {

  sim <- merge_simfiles(simno = i, indir = "data/", ftype = "min")
  sim <- truncate_sim(sim, at = 2600)
  vars.needed <- c(
    # HIV
    "incid", "hivtests.nprep", "ir100", "i.prev",

    # GC
    "incid.rgc", "incid.ugc",
    "incid.rgc.tttraj1", "incid.ugc.tttraj1",
    "incid.rgc.tttraj2", "incid.ugc.tttraj2",
    "ir100.gc", "ir100.gc.tttraj1", "ir100.gc.tttraj2",
    "prev.gc", "prev.gc.tttraj1", "prev.gc.tttraj2",
    "GCasympttests", "GCsympttests",
    "GCasympttests.tttraj1", "GCasympttests.tttraj2",
    "GCsympttests.tttraj1", "GCsympttests.tttraj2",
    "txGC", "txGC.tttraj1", "txGC.tttraj2",
    "txGC_asympt",
    "tt.traj.gc1", "tt.traj.gc2",

    # CT
    "incid.rct", "incid.uct",
    "incid.rct.tttraj1", "incid.uct.tttraj1",
    "incid.rct.tttraj2", "incid.uct.tttraj2",
    "ir100.ct", "ir100.ct.tttraj1", "ir100.ct.tttraj2",
    "prev.ct", "prev.ct.tttraj1", "prev.ct.tttraj2",
    "CTasympttests", "CTsympttests",
    "CTasympttests.tttraj1", "CTasympttests.tttraj2",
    "CTsympttests.tttraj1", "CTsympttests.tttraj2",
    "txCT", "txCT.tttraj1", "txCT.tttraj2",
    "txCT_asympt",
    "tt.traj.ct1", "tt.traj.ct2",

    # Combined
    "prev.sti", "prev.sti.tttraj1", "prev.sti.tttraj2",
    "tt.traj.sti1", "tt.traj.sti2",

    # Other
    "num")
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  filename <- paste0("data/sim.n", i, ".rda")
  save(sim, file = filename, compress = "gzip")

}
