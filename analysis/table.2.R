
## STI PrEP Table 2

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

names(sim$epi)
# "num.asympt.tx"         "num.asympt.cases"
# "num.asympt.tx.prep"    "num.asympt.cases.prep"
# "num.rect.tx"           "num.rect.cases"
# "num.rect.tx.prep"      "num.rect.cases.prep"

## base scenario
load("data/sim.n100.rda")
round(sim$param$prep.sti.screen.int * (12/52), 0)

inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx))
denom <- unname(colSums(sim$epi$num.asympt.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx))
denom <- unname(colSums(sim$epi$num.rect.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)


## Varying STI testing interval: n2000 to 2004

load("data/sim.n2004.rda")
round(sim$param$prep.sti.screen.int * (12/52), 0)

# total pop
inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx))
denom <- unname(colSums(sim$epi$num.asympt.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx))
denom <- unname(colSums(sim$epi$num.rect.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

# prep users
inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti.prep, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx.prep))
denom <- unname(colSums(sim$epi$num.asympt.cases.prep))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx.prep))
denom <- unname(colSums(sim$epi$num.rect.cases.prep))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)


## Varying prob STI treatment: n2005 to 2009
load("data/sim.n2008.rda")
round(sim$param$prep.sti.screen.int * (12/52), 0)
sim$param$prep.sti.prob.tx

# total pop
inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx))
denom <- unname(colSums(sim$epi$num.asympt.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx))
denom <- unname(colSums(sim$epi$num.rect.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

# prep users
inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti.prep, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx.prep))
denom <- unname(colSums(sim$epi$num.asympt.cases.prep))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx.prep))
denom <- unname(colSums(sim$epi$num.rect.cases.prep))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)



## Varying non-prep asymptomatic STI treatment: n2010 to 2014
load("data/sim.n2014.rda")
round(sim$param$prep.sti.screen.int * (12/52), 0)
sim$param$prep.sti.prob.tx
sim$param$gc.asympt.prob.tx

# total pop
inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx))
denom <- unname(colSums(sim$epi$num.asympt.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx))
denom <- unname(colSums(sim$epi$num.rect.cases))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

# prep users
inc <- as.numeric(colMeans(tail(sim$epi$ir100.sti.prep, 52)))
round(quantile(inc, c(0.5, 0.25, 0.75)), 2)

num <- unname(colSums(sim$epi$num.asympt.tx.prep))
denom <- unname(colSums(sim$epi$num.asympt.cases.prep))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)

num <- unname(colSums(sim$epi$num.rect.tx.prep))
denom <- unname(colSums(sim$epi$num.rect.cases.prep))
vec <- num / denom
round(quantile(vec, c(0.5, 0.25, 0.75)), 3)




