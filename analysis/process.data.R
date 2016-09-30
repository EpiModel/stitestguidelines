
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
fn <- list.files("data/", pattern = "n10[0-9][0-9]", full.names = TRUE)[-1]
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  save(sim, file = i, compress = TRUE)
  cat("*")
}

# Remove unneeded variables n1000 data
fn <- list.files("data/", pattern = "n10[0-9][0-9]", full.names = TRUE)
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
fn <- list.files("data/", pattern = "n201[0-9]", full.names = TRUE)
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  save(sim, file = i, compress = TRUE)
  cat("*")
}
