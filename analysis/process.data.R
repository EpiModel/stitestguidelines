
## Process STI PrEP Data

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# unlink("data/*")
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")
system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n200[0-9].rda data/")
system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n201[0-9].rda data/")

( fn <- list.files("data/", full.names = TRUE) )

# Truncate n1000 data
fn <- list.files("data/", pattern = "n10[0-9][0-9]", full.names = TRUE)
fn <- list.files("data/", pattern = "n100.rda", full.names = TRUE)
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  save(sim, file = i, compress = TRUE)
  cat("*")
}

# Remove unneeded variables n1000 data
fn <- list.files("data/", pattern = "n10[0-9][0-9]", full.names = TRUE)
fn <- list.files("data/", pattern = "n100.rda", full.names = TRUE)
for (i in fn) {
  load(i)
  vars.needed <- c("ir100", "incid", "ir100.gc", "incid.gc", "ir100.ct", "incid.ct",
                   "prepCurr", "i.prev", "prev.gc", "prev.ct")
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  save(sim, file = i, compress = "xz")
  cat("\n sim ", i, " complete")
}

# Truncate n2000 data
fn <- list.files("data/", pattern = "n200[0-9]", full.names = TRUE)
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  save(sim, file = i, compress = TRUE)
  cat("*")
}

# Remove unneeded variables n2000 data
fn <- list.files("data/", pattern = "n200[0-9]", full.names = TRUE)
for (i in fn) {
  load(i)
  vars.needed <- c("ir100", "incid", "ir100.gc", "incid.gc", "ir100.ct", "incid.ct",
                   "prepCurr", "i.prev", "prev.gc", "prev.ct",
                   "ir100.sti", "ir100.sti.prep", "num.asympt.tx",
                   "num.asympt.cases", "num.rect.tx", "num.rect.cases")
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  save(sim, file = i, compress = "xz")
  cat("\n", i, "complete")
}


# truncate and limit sim 3000/4000 files
suppressMessages(library("EpiModelHIV"))
fn <- list.files(pattern = "n[3-4][0-9][0-9][0-9].rda")
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  vars.needed <- c("ir100", "incid", "ir100.gc", "incid.gc",
                   "ir100.ct", "incid.ct", "prepCurr",
                   "ir100.sti", "ir100.sti.prep")
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  out.fn <- paste0("save/", i)
  save(sim, file = out.fn, compress = "xz")
  file.remove(i)
  cat(i, "\n")
}
