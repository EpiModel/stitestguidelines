
## STI PrEP Table 3: HIV outcomes


rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")


load("data/sim.n100.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying PrEP coverage 10%, 20%, 40%: 1036, 1037, 1039
load("data/sim.n1039.rda")
sim$param$prep.coverage
sim$param$rcomp.prob
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)


# Vary risk compensation 0%, 40%, 100%: 1003, 1039, 1093
load("data/sim.n1093.rda")
sim$param$prep.coverage
sim$param$rcomp.prob
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Vary STI testing interval 1, 3, 6 months: 2000, 2001, 1039
load("data/sim.n2001.rda")
sim$param$prep.coverage
sim$param$rcomp.prob
round(sim$param$prep.sti.screen.int * (12/52), 0)
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)


