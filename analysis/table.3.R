## STI PrEP Table 2
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
load("data/followup/sim.n3000.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying Lower-Risk
# 3131:3141, 3054, Lower-Risk = 40%, 28 days to 364 days by ~28 days, Higher-Risk = 0%, 182 days
# 28, 63, 91, 119, 147, 182, 210, 238, 273, 301, 329, 364

load("data/followup/sim.n3131.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3132.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3133.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3134.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3135.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3136.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3137.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3138.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3139.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3140.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3141.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3054.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying Higher-Risk
# 3153, 3156, 3158, 3160, 3162, 3164, 3166, 3168, 3170, 3172, 3014, Annual = 0%, 364 days, Higher-Risk = 40%, 182 days
# 28, 49, 63, 77, 91, 105, 119, 133, 147, 161, 182 

load("data/followup/sim.n3153.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3156.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3158.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3160.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3162.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3164.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3166.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3168.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3170.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3172.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3014.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)


