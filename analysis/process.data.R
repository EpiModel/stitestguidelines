## Process STI Testing Guidelines Data

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# unlink("data/*")
# system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")
# system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n200[0-9].rda data/")
# system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n201[0-9].rda data/")

( fn <- list.files("data/followup", full.names = TRUE) )
# fn <- list.files(pattern = "data/followup/n[3-4][0-9][0-9][0-9].rda")
for (i in fn) {
    load(i)
    sim <- truncate_sim(sim, at = 2600)
    save(sim, file = i, compress = TRUE)
    cat("*")
}

# truncate and limit sim 3000/4000 files on Hyak
suppressMessages(library("EpiModelHIV"))
fn <- list.files(pattern = "n[3-4][0-9][0-9][0-9].rda")
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  vars.needed <- c("ir100", "incid", "ir100.gc", "incid.gc",
                   "ir100.ct", "incid.ct", "ir100.syph", "incid.syph",
                   "ir100.sti", "ir100.sti.prep", "totalstiasympttests",
                   "totalGCasympttests", "totalCTasympttests", "totalsyphasympttests",
                   "i.prev", "prev.primsecosyph", "prev.syph", 
                   "totalhivtests", "totalhivtests.prep","time.on.prep", 
                   "stage.time.ar.ndx", "stage.time.af.ndx", "stage.time.chronic.ndx", "stage.time.aids.ndx",
                   "stage.time.ar.dx", "stage.time.af.dx", "stage.time.chronic.dx", "stage.time.aids.dx",
                   "stage.time.ar.art", "stage.time.af.art", "stage.time.chronic.art", "stage.time.aids.art",
                   "totalrGCasympttests.prep", "totaluGCasympttests.prep", "totalGCasympttests.prep",
                   "totalrCTasympttests.prep", "totaluCTasympttests.prep", "totalCTasympttests.prep",
                   "totalsyphasympttests.prep")
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  out.fn <- paste0("data/followup/", i)
  save(sim, file = out.fn, compress = "gzip")
  file.remove(i)
  cat(i, "\n")
}
