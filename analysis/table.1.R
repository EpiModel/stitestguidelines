## STI PrEP Table 1
# Varying Indications for High-Risk Testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No testing
load("data/sim.n3000.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

## Varying Indications:

# 3142 - STI
load("data/sim.n3142.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3143 - recent partners
load("data/sim.n3143.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3144 - new partners
load("data/sim.n3144.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3145 - partner who has multiple partners
load("data/sim.n3145.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3146 - partner with a STI
load("data/sim.n3146.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3147 - any CAI in a non-main
load("data/sim.n3147.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3148 - any CAI
load("data/sim.n3148.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3149 - recent or new partners
load("data/sim.n3149.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3150 - sti, recent, or new partners
load("data/sim.n3150.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3151 - CAI in non-main or any CAI
load("data/sim.n3151.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3152 - partner with multiple partners or with a STI
load("data/sim.n3152.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)