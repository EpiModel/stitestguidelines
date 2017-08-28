## Process STI Testing Guidelines Data

# Bulk on Hyak -----------------------------------------------------------------
rm(list = ls())
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
fn <- list.files(pattern = "n[3-5][0-9][0-9][0-9].rda")
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 5200)
  vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                   "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                   "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                   "ir100.sti", "ir100.sti.prep",
                   "GCasympttests", "uGCasympttests", "rGCasympttests",
                   "CTasympttests", "uCTasympttests", "rCTasympttests",
                   "syphasympttests", "stiasympttests",
                   "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos",
                   "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos",
                   "syphasympttests.pos", "syphearlyasympttests.pos",
                   "syphlateasympttests.pos", "stiasympttests.pos",
                   "GCasympttests.hivneg", "uGCasympttests.hivneg", "rGCasympttests.hivneg",
                   "CTasympttests.hivneg", "uCTasympttests.hivneg", "rCTasympttests.hivneg",
                   "syphasympttests.hivneg", "stiasympttests.hivneg",
                   "GCasympttests.pos.hivneg", "uGCasympttests.pos.hivneg", "rGCasympttests.pos.hivneg",
                   "CTasympttests.pos.hivneg", "uCTasympttests.pos.hivneg", "rCTasympttests.pos.hivneg",
                   "syphasympttests.pos.hivneg", "syphearlyasympttests.pos.hivneg",
                   "syphlateasympttests.pos.hivneg", "stiasympttests.pos.hivneg",
                   "GCasympttests.hivpos", "uGCasympttests.hivpos", "rGCasympttests.hivpos",
                   "CTasympttests.hivpos", "uCTasympttests.hivpos", "rCTasympttests.hivpos",
                   "syphasympttests.hivpos", "stiasympttests.hivpos",
                   "GCasympttests.pos.hivpos", "uGCasympttests.pos.hivpos", "rGCasympttests.pos.hivpos",
                   "CTasympttests.pos.hivpos", "uCTasympttests.pos.hivpos", "rCTasympttests.pos.hivpos",
                   "syphasympttests.pos.hivpos", "syphearlyasympttests.pos.hivpos",
                   "syphlateasympttests.pos.hivpos", "stiasympttests.pos.hivpos",
                   "GCsympttests", "uGCsympttests", "rGCsympttests",
                   "CTsympttests", "uCTsympttests", "rCTsympttests",
                   "syphsympttests", "stisympttests",
                   "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                   "i.prev", "prev.primsecosyph", "prev.syph",
                   "prev.gc", "prev.rgc", "prev.ugc",
                   "prev.ct", "prev.rct", "prev.uct",
                   "prev.rgc.hivneg.only", "prev.ugc.hivneg.only",
                   "prev.gc.hivneg.only", "prev.rct.hivneg.only",
                   "prev.uct.hivneg.only", "prev.ct.hivneg.only",
                   "prev.syph.hivneg.only", "prev.rgc.hivpos.only",
                   "prev.ugc.hivpos.only", "prev.gc.hivpos.only",
                   "prev.rct.hivpos.only", "prev.uct.hivpos.only",
                   "prev.ct.hivpos.only", "prev.syph.hivpos.only",
                   "prev.hivposmultsti", "prev.hivnegmultsti",
                   "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT", "txasympt",
                   "txGC_asympt", "txCT_asympt", "txsyph_asympt",
                   "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                   "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                   "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                   "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                   "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                   "deathage", "stiactiveind.prop", "stiactiveind",
                   "recentpartners", "recentpartners.prop",
                   "num.earlydiagsyph", "num.latediagsyph", "num.newearlydiagsyph",
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
  #file.remove(i)
  cat(i, "\n")
}


### 1 by 1 processing on Hyak - cd /gscratch/csde/kweiss2/sti/data -------------
rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
load("sim.n4009.rda")
sim <- truncate_sim(sim, at = 5200)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.sti.prep",
                 "GCasympttests", "uGCasympttests", "rGCasympttests",
                 "CTasympttests", "uCTasympttests", "rCTasympttests",
                 "syphasympttests", "stiasympttests",
                 "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos",
                 "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos",
                 "syphasympttests.pos", "syphearlyasympttests.pos",
                 "syphlateasympttests.pos", "stiasympttests.pos",
                 "GCasympttests.hivneg", "uGCasympttests.hivneg", "rGCasympttests.hivneg",
                 "CTasympttests.hivneg", "uCTasympttests.hivneg", "rCTasympttests.hivneg",
                 "syphasympttests.hivneg", "stiasympttests.hivneg",
                 "GCasympttests.pos.hivneg", "uGCasympttests.pos.hivneg", "rGCasympttests.pos.hivneg",
                 "CTasympttests.pos.hivneg", "uCTasympttests.pos.hivneg", "rCTasympttests.pos.hivneg",
                 "syphasympttests.pos.hivneg", "syphearlyasympttests.pos.hivneg",
                 "syphlateasympttests.pos.hivneg", "stiasympttests.pos.hivneg",
                 "GCasympttests.hivpos", "uGCasympttests.hivpos", "rGCasympttests.hivpos",
                 "CTasympttests.hivpos", "uCTasympttests.hivpos", "rCTasympttests.hivpos",
                 "syphasympttests.hivpos", "stiasympttests.hivpos",
                 "GCasympttests.pos.hivpos", "uGCasympttests.pos.hivpos", "rGCasympttests.pos.hivpos",
                 "CTasympttests.pos.hivpos", "uCTasympttests.pos.hivpos", "rCTasympttests.pos.hivpos",
                 "syphasympttests.pos.hivpos", "syphearlyasympttests.pos.hivpos",
                 "syphlateasympttests.pos.hivpos", "stiasympttests.pos.hivpos",
                 "GCsympttests", "uGCsympttests", "rGCsympttests",
                 "CTsympttests", "uCTsympttests", "rCTsympttests",
                 "syphsympttests", "stisympttests",
                 "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                 "i.prev", "prev.primsecosyph", "prev.syph",
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct",
                 "prev.rgc.hivneg.only", "prev.ugc.hivneg.only",
                 "prev.gc.hivneg.only", "prev.rct.hivneg.only",
                 "prev.uct.hivneg.only", "prev.ct.hivneg.only",
                 "prev.syph.hivneg.only", "prev.rgc.hivpos.only",
                 "prev.ugc.hivpos.only", "prev.gc.hivpos.only",
                 "prev.rct.hivpos.only", "prev.uct.hivpos.only",
                 "prev.ct.hivpos.only", "prev.syph.hivpos.only",
                 "prev.hivposmultsti", "prev.hivnegmultsti",
                 "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT", "txasympt",
                 "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                 "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                 "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                 "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                 "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                 "deathage", "stiactiveind.prop", "stiactiveind",
                 "recentpartners", "recentpartners.prop",
                 "num.earlydiagsyph", "num.latediagsyph", "num.newearlydiagsyph",
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
save(sim, file = "followup/sim.n4009.rda", compress = "gzip")


## Locally merge files --------------------------------------------------------
sim <- merge_simfiles(3118, indir = "data/", ftype = "min")
sim <- truncate_sim(sim, at = 5200)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.sti.prep",
                 "GCasympttests", "uGCasympttests", "rGCasympttests",
                 "CTasympttests", "uCTasympttests", "rCTasympttests",
                 "syphasympttests", "stiasympttests",
                 "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos",
                 "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos",
                 "syphasympttests.pos", "syphearlyasympttests.pos",
                 "syphlateasympttests.pos", "stiasympttests.pos",
                 "GCasympttests.hivneg", "uGCasympttests.hivneg", "rGCasympttests.hivneg",
                 "CTasympttests.hivneg", "uCTasympttests.hivneg", "rCTasympttests.hivneg",
                 "syphasympttests.hivneg", "stiasympttests.hivneg",
                 "GCasympttests.pos.hivneg", "uGCasympttests.pos.hivneg", "rGCasympttests.pos.hivneg",
                 "CTasympttests.pos.hivneg", "uCTasympttests.pos.hivneg", "rCTasympttests.pos.hivneg",
                 "syphasympttests.pos.hivneg", "syphearlyasympttests.pos.hivneg",
                 "syphlateasympttests.pos.hivneg", "stiasympttests.pos.hivneg",
                 "GCasympttests.hivpos", "uGCasympttests.hivpos", "rGCasympttests.hivpos",
                 "CTasympttests.hivpos", "uCTasympttests.hivpos", "rCTasympttests.hivpos",
                 "syphasympttests.hivpos", "stiasympttests.hivpos",
                 "GCasympttests.pos.hivpos", "uGCasympttests.pos.hivpos", "rGCasympttests.pos.hivpos",
                 "CTasympttests.pos.hivpos", "uCTasympttests.pos.hivpos", "rCTasympttests.pos.hivpos",
                 "syphasympttests.pos.hivpos", "syphearlyasympttests.pos.hivpos",
                 "syphlateasympttests.pos.hivpos", "stiasympttests.pos.hivpos",
                 "GCsympttests", "uGCsympttests", "rGCsympttests",
                 "CTsympttests", "uCTsympttests", "rCTsympttests",
                 "syphsympttests", "stisympttests",
                 "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                 "i.prev", "prev.primsecosyph", "prev.syph",
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct",
                 "prev.rgc.hivneg.only", "prev.ugc.hivneg.only",
                 "prev.gc.hivneg.only", "prev.rct.hivneg.only",
                 "prev.uct.hivneg.only", "prev.ct.hivneg.only",
                 "prev.syph.hivneg.only", "prev.rgc.hivpos.only",
                 "prev.ugc.hivpos.only", "prev.gc.hivpos.only",
                 "prev.rct.hivpos.only", "prev.uct.hivpos.only",
                 "prev.ct.hivpos.only", "prev.syph.hivpos.only",
                 "prev.hivposmultsti", "prev.hivnegmultsti",
                 "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT", "txasympt",
                 "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                 "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                 "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                 "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                 "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                 "deathage", "stiactiveind.prop", "stiactiveind",
                 "recentpartners", "recentpartners.prop",
                 "num.earlydiagsyph", "num.latediagsyph", "num.newearlydiagsyph",
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

#### Merge on Hyak ------------------------------------------------------------
sim <- merge_simfiles(3003, ftype = "min")
sim <- truncate_sim(sim, at = 2600)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.sti.prep",
                 "GCasympttests", "uGCasympttests", "rGCasympttests",
                 "CTasympttests", "uCTasympttests", "rCTasympttests",
                 "syphasympttests", "stiasympttests",
                 "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos",
                 "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos",
                 "syphasympttests.pos", "syphearlyasympttests.pos",
                 "syphlateasympttests.pos", "stiasympttests.pos",
                 "GCasympttests.hivneg", "uGCasympttests.hivneg", "rGCasympttests.hivneg",
                 "CTasympttests.hivneg", "uCTasympttests.hivneg", "rCTasympttests.hivneg",
                 "syphasympttests.hivneg", "stiasympttests.hivneg",
                 "GCasympttests.pos.hivneg", "uGCasympttests.pos.hivneg", "rGCasympttests.pos.hivneg",
                 "CTasympttests.pos.hivneg", "uCTasympttests.pos.hivneg", "rCTasympttests.pos.hivneg",
                 "syphasympttests.pos.hivneg", "syphearlyasympttests.pos.hivneg",
                 "syphlateasympttests.pos.hivneg", "stiasympttests.pos.hivneg",
                 "GCasympttests.hivpos", "uGCasympttests.hivpos", "rGCasympttests.hivpos",
                 "CTasympttests.hivpos", "uCTasympttests.hivpos", "rCTasympttests.hivpos",
                 "syphasympttests.hivpos", "stiasympttests.hivpos",
                 "GCasympttests.pos.hivpos", "uGCasympttests.pos.hivpos", "rGCasympttests.pos.hivpos",
                 "CTasympttests.pos.hivpos", "uCTasympttests.pos.hivpos", "rCTasympttests.pos.hivpos",
                 "syphasympttests.pos.hivpos", "syphearlyasympttests.pos.hivpos",
                 "syphlateasympttests.pos.hivpos", "stiasympttests.pos.hivpos",
                 "GCsympttests", "uGCsympttests", "rGCsympttests",
                 "CTsympttests", "uCTsympttests", "rCTsympttests",
                 "syphsympttests", "stisympttests",
                 "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                 "i.prev", "prev.primsecosyph", "prev.syph",
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct",
                 "prev.rgc.hivneg.only", "prev.ugc.hivneg.only",
                 "prev.gc.hivneg.only", "prev.rct.hivneg.only",
                 "prev.uct.hivneg.only", "prev.ct.hivneg.only",
                 "prev.syph.hivneg.only", "prev.rgc.hivpos.only",
                 "prev.ugc.hivpos.only", "prev.gc.hivpos.only",
                 "prev.rct.hivpos.only", "prev.uct.hivpos.only",
                 "prev.ct.hivpos.only", "prev.syph.hivpos.only",
                 "prev.hivposmultsti", "prev.hivnegmultsti",
                 "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT", "txasympt",
                 "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                 "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                 "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                 "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                 "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                 "deathage", "stiactiveind.prop", "stiactiveind",
                 "recentpartners", "recentpartners.prop",
                 "num.earlydiagsyph", "num.latediagsyph", "num.newearlydiagsyph",
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
save(sim, file = "followup/sim.3003.rda", compress = "gzip")
