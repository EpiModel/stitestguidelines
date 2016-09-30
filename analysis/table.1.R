
## STI PrEP Table 1

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")


load("data/sim.n100.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying PrEP coverage 10% to 90% at 40% risk comp: n1036 to 1044
load("data/sim.n1044.rda")
sim$param$prep.coverage
sim$param$rcomp.prob
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)


# Varying Risk Comp 0% to 100% at 40% coverage: 1003, 1012, 1021, 1030, 1039,
#                                               1048, 1057, 1066, 1075, 1084, 1093
load("data/sim.n1093.rda")
sim$param$prep.coverage
sim$param$rcomp.prob
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)
