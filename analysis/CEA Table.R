## STI PrEP CEA Table
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
# load("data/followup/sim.n3000.rda")
# sim.base <- sim
# epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Newer way
# Base - No annual or high-risk
# haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
# ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
# incid.base <- unname(colSums(sim.base$epi$incid))
# 
# haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
# ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
# incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
# 
# haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
# ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
# incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
# 
# haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
# ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
# incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))

## Base STI lower-risk testing interval (364 days): n3054
## Varying STI lower-risk testing interval: n3131-n3141
## Base STI higher-risk testing interval: n3014
## Varying STI higher-risk testing interval: n3153 - n3173
#sims <- c(3000, 3131:3141, 3054, 3153, 3156, 3158, 3160, 3162, 3164, 3166, 3168, 3170, 3172, 3014)

sims <- c(3000, 3131:3141)

qnt.low <- 0.25
qnt.high <- 0.75

annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
QALY <- rep(NA, length(sims))
hivtestcosts <- rep(NA, length(sims))
hivhealthcosts <- rep(NA, length(sims))
stisympttestcosts <- rep(NA, length(sims))
gc.asympttestcosts <- rep(NA, length(sims))
ct.asympttestcosts <- rep(NA, length(sims))
syph.asympttestcosts <- rep(NA, length(sims))
rect.asympttestcosts <- rep(NA, length(sims))
ureth.asympttestcosts <- rep(NA, length(sims))
stiasympttestcosts <- rep(NA, length(sims))
stitxcosts <- rep(NA, length(sims))

# add sims to data frame as an object?
df <- data.frame(annint, hrint, anncov, hrcov, QALY, hivtestcosts, hivhealthcosts, stisympttestcosts,
                 gc.asympttestcosts, ct.asympttestcosts, syph.asympttestcosts,
                 rect.asympttestcosts, ureth.asympttestcosts, stitxcosts,
                 stiasympttestcosts)

for (i in seq_along(sims)) {
    
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    
    df$annint[i] <- sim$param$stitest.active.int
    df$hrint[i] <- sim$param$sti.highrisktest.int
    df$anncov[i] <- sim$param$stianntest.coverage
    df$hrcov[i] <- sim$param$stihighrisktest.coverage
    
    # Time in HIV stages
    time.hivneg <- unname(colSums(sim$epi$time.hivneg) / 52)
    time.on.prep <- unname(colSums(sim$epi$time.on.prep) / 52)
    time.off.prep <- unname(colSums(sim$epi$time.off.prep) / 52)
    stage.time.ar.ndx <- unname(colSums(sim$epi$stage.time.ar.ndx) / 52)
    stage.time.af.ndx <- unname(colSums(sim$epi$stage.time.af.ndx) / 52)
    stage.time.chronic.ndx <- unname(colSums(sim$epi$stage.time.chronic.ndx) / 52)
    stage.time.aids.ndx <- unname(colSums(sim$epi$stage.time.aids.ndx) / 52)
    stage.time.ar.dx <- unname(colSums(sim$epi$stage.time.ar.dx) / 52)
    stage.time.af.dx <- unname(colSums(sim$epi$stage.time.af.dx) / 52)
    stage.time.chronic.dx <- unname(colSums(sim$epi$stage.time.chronic.dx) / 52)
    stage.time.aids.dx <- unname(colSums(sim$epi$stage.time.aids.dx) / 52)
    stage.time.ar.art <- unname(colSums(sim$epi$stage.time.ar.art) / 52)
    stage.time.af.art <- unname(colSums(sim$epi$stage.time.af.art) / 52)
    stage.time.chronic.art <- unname(colSums(sim$epi$stage.time.chronic.art) / 52)
    stage.time.aids.art <- unname(colSums(sim$epi$stage.time.aids.art) / 52)
    
    # Total tests
    totalhivtests <- unname(sim$epi$totalhivtests[520, ])
    totalrGCsympttests <- unname(sim$epi$totalrGCsympttests[520, ])
    totaluGCsympttests <- unname(sim$epi$totaluGCsympttests[520, ])
    totalGCsympttests <- unname(sim$epi$totalGCsympttests[520, ])
    totalrCTsympttests <- unname(sim$epi$totalrCTsympttests[520, ])
    totaluCTsympttests <- unname(sim$epi$totaluCTsympttests[520, ])
    totalCTsympttests <- unname(sim$epi$totalCTsympttests[520, ])
    totalsyphsympttests <- unname(sim$epi$totalsyphsympttests[520, ])
    totalrGCasympttests <- unname(sim$epi$totalrGCasympttests[520, ])
    totaluGCasympttests <- unname(sim$epi$totaluGCasympttests[520, ])
    totalGCasympttests <- unname(sim$epi$totalGCasympttests[520, ])
    totalrCTasympttests <- unname(sim$epi$totalrCTasympttests[520, ])
    totaluCTasympttests <- unname(sim$epi$totaluCTasympttests[520, ])
    totalCTasympttests <- unname(sim$epi$totalCTasympttests[520, ])
    totalsyphasympttests <- unname(sim$epi$totalsyphasympttests[520, ])
    totalhivtests.pos <- unname(sim$epi$totalhivtests.pos[520, ])
    totalrGCasympttests.pos <- unname(sim$epi$totalrGCasympttests.pos[520, ])
    totaluGCasympttests.pos <- unname(sim$epi$totaluGCasympttests.pos[520, ])
    totalGCasympttests.pos <- unname(sim$epi$totalGCasympttests.pos[520, ])
    totalrCTasympttests.pos <- unname(sim$epi$totalrCTasympttests.pos[520, ])
    totaluCTasympttests.pos <- unname(sim$epi$totaluCTasympttests.pos[520, ])
    totalCTasympttests.pos <- unname(sim$epi$totalCTasympttests.pos[520, ])
    totalsyphasympttests.pos <- unname(sim$epi$totalsyphasympttests.pos[520, ])
    
    # STI Tx
    txCT <- unname(colSums(sim$epi$txCT))
    txGC <- unname(colSums(sim$epi$txGC))
    txsyph <- unname(colSums(sim$epi$txsyph))
    txearlysyph <- unname(colSums(sim$epi$txearlysyph))
    txlatesyph <- unname(colSums(sim$epi$txlatesyph))
    
    # QALY
    df$QALY[i] <- sum((time.hivneg * 1), 
                      (stage.time.ar.ndx * 0.865), (stage.time.af.ndx * 0.865),
                      (stage.time.chronic.ndx * 0.865), (stage.time.aids.ndx * 0.72),
                      (stage.time.ar.dx * 0.795), (stage.time.af.dx * 0.795),
                      (stage.time.chronic.dx * 0.795), (stage.time.aids.dx * 0.72),
                      (stage.time.ar.art * 0.83), (stage.time.af.art * 0.83),
                      (stage.time.chronic.art * 0.83), (stage.time.aids.art * 0.82))

    df$hivtestcosts[i] <- sum(((totalhivtests - totalhivtests.pos) * 64.75), (totalhivtests.pos * 614.75))
    df$healthcosts[i] <- sum((time.hivneg * 4469.81), 
                             (stage.time.ar.ndx * 4502.83), (stage.time.ar.dx * 4502.83), (stage.time.af.art * 17071.8),
                             (stage.time.af.ndx * 4502.83),  (stage.time.af.dx * 4502.83), (stage.time.af.art * 17071.8),
                             (stage.time.chronic.ndx * 10558.70), (stage.time.chronic.dx * 10558.70), (stage.time.chronic.art * 23842),
                             (stage.time.aids.ndx * 28533.69), (stage.time.aids.dx * 28533.69), (stage.time.aids.art * 27990.41))
    df$stisympttestcosts[i] <- sum((totalrGCsympttests * 45.62), (totaluGCsympttests * 45.62),
                                   (totalrCTsympttests * 45.62), (totaluCTsympttests * 45.62),
                                   (totalsyphsympttests * 33.96))
    df$gc.asympttestcosts[i] <- sum((totalGCasympttests * 45.62))
    df$ct.asympttestcosts[i] <- sum((totalCTasympttests * 45.62))
    df$syph.asympttestcosts[i] <- sum((totalsyphasympttests * 33.96))
    df$rect.asympttestcosts[i] <- sum((totalrCTasympttests * 45.62), (totalrGCasympttests * 45.62))
    df$ureth.asympttestcosts[i] <- sum((totaluCTasympttests * 45.62), (totaluGCasympttests * 45.62))
    df$stiasympttestcosts[i] <- sum(sum((totalrGCasympttests * 45.62), (totaluGCasympttests * 45.62),
                                        (totalrCTasympttests * 45.62), (totaluCTasympttests * 45.62),
                                        (totalsyphasympttests * 33.96)))
    df$stitxcosts[i] <- sum((txlatesyph * 99.35 * 3), (txearlysyph * 99.35), (txGC * 53.35), (txCT * 53.35))

    cat("*")
    
}

View(df)
