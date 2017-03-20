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
# 3014, 3025, 3036, 3047, 3058, 3069, 3080, 3091, 3102, 3113, 3124: Annual = 0.1 - 1.0 by 0.1, 364 days, HR = 40%, 182 days
load("data/followup/sim.n3014.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3025.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3036.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3047.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3058.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3069.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3080.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3091.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3102.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3113.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3124.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying Higher-Risk
# 3054 - 3064 Annual = 40%, 364 days, HR = 0.0 - 1.0 by 0.1, 182 days
load("data/followup/sim.n3054.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3055.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3056.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3057.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

## Same as above scenario
load("data/followup/sim.n3058.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3059.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3060.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3061.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3062.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3063.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3064.rda")
sim$param$stianntest.coverage
sim$param$stihrtest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

