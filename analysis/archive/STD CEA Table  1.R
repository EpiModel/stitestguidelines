## STI Testing Indications CEA Table 1
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Varying Lower-Risk Coverage
# 3012, 3023, 3034, 3045, 3056, 3067, 3078, 3089, 3100, 3111, 3122: Annual = 0.0 - 1.0 by 0.1, 364 days, HR = 20%, 182 days
# Varying Higher-Risk Coverage
# 3032:3042: Higher-risk = 0.0 - 1.0 by 0.1, 182 days, Ann = 20%, 364 days
sims <- c(3000, 3012, 3023, 3034, 3045, 3056, 3067, 3078, 3089, 3100, 3111, 3122, 3032:3042)

qnt.low <- 0.25
qnt.high <- 0.75

partcut <- rep(NA, length(sims))
probtx <- rep(NA, length(sims))
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
stiasympttests <- rep(NA, length(sims))
stitxcosts <- rep(NA, length(sims))
overallcosts <- rep(NA, length(sims))
incremover <- rep(NA, length(sims))
incremQALYover <- rep(NA, length(sims))
costpqaly <- rep(NA, length(sims))
costpoverqaly <- rep(NA, length(sims))

# add sims to data frame as an object?
df <- data.frame(partcutprobtx, annint, hrint, anncov, hrcov, QALY, stiasympttests, hivtestcosts, hivhealthcosts, stisympttestcosts,
                 gc.asympttestcosts, ct.asympttestcosts, syph.asympttestcosts,
                 rect.asympttestcosts, ureth.asympttestcosts, stitxcosts,
                 stiasympttestcosts, overallcosts, incremover,
                 incremQALYover, costpqaly, costpoverqaly)

for (i in seq_along(sims)) {

    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)

    #sim <- truncate_sim(sim, at = 2600)
    df$partcut[i] <- sim$param$partnercutoff
    df$probtx[i] <- sim$param$gc.asympt.prob.tx
    df$annint[i] <- sim$param$stitest.active.int
    df$hrint[i] <- sim$param$sti.highrisktest.int
    df$anncov[i] <- sim$param$stianntest.coverage
    df$hrcov[i] <- sim$param$stihighrisktest.coverage

    # Time in HIV stages
    time.hivneg <- ((quantile(unname(colSums(sim$epi$time.hivneg[2:53,])), probs = 0.5) * (0.97^0)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[54:105,])), probs = 0.5) * (0.97^1)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[106:157,])), probs = 0.5) * (0.97^2)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[158:209,])), probs = 0.5) * (0.97^3)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[210:261,])), probs = 0.5) * (0.97^4)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[262:313,])), probs = 0.5) * (0.97^5)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[314:365,])), probs = 0.5) * (0.97^6)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[366:417,])), probs = 0.5) * (0.97^7)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[418:469,])), probs = 0.5) * (0.97^8)) +
                      (quantile(unname(colSums(sim$epi$time.hivneg[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.ar.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.ar.ndx[2:53,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[54:105,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[106:157,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[158:209,])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[210:261,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[262:313,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[314:365,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[366:417,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[418:469,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.ar.dx <- ((quantile(unname(colSums(sim$epi$stage.time.ar.dx[2:53,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[54:105,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[106:157,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[158:209,])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[210:261,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[262:313,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[314:365,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[366:417,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[418:469,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.ar.dx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.af.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.af.ndx[2:53,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[54:105,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[106:157,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[158:209,])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[210:261,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[262:313,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[314:365,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[366:417,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[418:469,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.af.ndx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.af.dx <- ((quantile(unname(colSums(sim$epi$stage.time.af.dx[2:53,])), probs = 0.5) * (0.97^0)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[54:105,])), probs = 0.5) * (0.97^1)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[106:157,])), probs = 0.5) * (0.97^2)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[158:209,])), probs = 0.5) * (0.97^3)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[210:261,])), probs = 0.5) * (0.97^4)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[262:313,])), probs = 0.5) * (0.97^5)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[314:365,])), probs = 0.5) * (0.97^6)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[366:417,])), probs = 0.5) * (0.97^7)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[418:469,])), probs = 0.5) * (0.97^8)) +
                           (quantile(unname(colSums(sim$epi$stage.time.af.dx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.early.chronic.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[2:53,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[54:105,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[106:157,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[158:209,])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[210:261,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[262:313,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[314:365,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[366:417,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[418:469,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.early.chronic.dx.yrone <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[2:53,])), probs = 0.5) * (0.97^0)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[54:105,])), probs = 0.5) * (0.97^1)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[106:157,])), probs = 0.5) * (0.97^2)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[158:209,])), probs = 0.5) * (0.97^3)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[210:261,])), probs = 0.5) * (0.97^4)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[262:313,])), probs = 0.5) * (0.97^5)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[314:365,])), probs = 0.5) * (0.97^6)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[366:417,])), probs = 0.5) * (0.97^7)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[418:469,])), probs = 0.5) * (0.97^8)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.early.chronic.dx.yrstwotolate <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[2:53,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[54:105,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[106:157,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[158:209,])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[210:261,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[262:313,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[314:365,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[366:417,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[418:469,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.early.chronic.art <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[2:53,])), probs = 0.5) * (0.97^0)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[54:105,])), probs = 0.5) * (0.97^1)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[106:157,])), probs = 0.5) * (0.97^2)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[158:209,])), probs = 0.5) * (0.97^3)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[210:261,])), probs = 0.5) * (0.97^4)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[262:313,])), probs = 0.5) * (0.97^5)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[314:365,])), probs = 0.5) * (0.97^6)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[366:417,])), probs = 0.5) * (0.97^7)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[418:469,])), probs = 0.5) * (0.97^8)) +
                           (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.late.chronic.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[2:53,])), probs = 0.5) * (0.97^0)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[54:105,])), probs = 0.5) * (0.97^1)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[106:157,])), probs = 0.5) * (0.97^2)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[158:209,])), probs = 0.5) * (0.97^3)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[210:261,])), probs = 0.5) * (0.97^4)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[262:313,])), probs = 0.5) * (0.97^5)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[314:365,])), probs = 0.5) * (0.97^6)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[366:417,])), probs = 0.5) * (0.97^7)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[418:469,])), probs = 0.5) * (0.97^8)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.late.chronic.dx <- ((quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[2:53,])), probs = 0.5) * (0.97^0)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[54:105,])), probs = 0.5) * (0.97^1)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[106:157,])), probs = 0.5) * (0.97^2)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[158:209,])), probs = 0.5) * (0.97^3)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[210:261,])), probs = 0.5) * (0.97^4)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[262:313,])), probs = 0.5) * (0.97^5)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[314:365,])), probs = 0.5) * (0.97^6)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[366:417,])), probs = 0.5) * (0.97^7)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[418:469,])), probs = 0.5) * (0.97^8)) +
                                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.late.chronic.art <- ((quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[2:53,])), probs = 0.5) * (0.97^0)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[54:105,])), probs = 0.5) * (0.97^1)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[106:157,])), probs = 0.5) * (0.97^2)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[158:209,])), probs = 0.5) * (0.97^3)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[210:261,])), probs = 0.5) * (0.97^4)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[262:313,])), probs = 0.5) * (0.97^5)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[314:365,])), probs = 0.5) * (0.97^6)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[366:417,])), probs = 0.5) * (0.97^7)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[418:469,])), probs = 0.5) * (0.97^8)) +
                                       (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.aids.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.aids.ndx[2:53,])), probs = 0.5) * (0.97^0)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[54:105,])), probs = 0.5) * (0.97^1)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[106:157,])), probs = 0.5) * (0.97^2)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[158:209,])), probs = 0.5) * (0.97^3)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[210:261,])), probs = 0.5) * (0.97^4)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[262:313,])), probs = 0.5) * (0.97^5)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[314:365,])), probs = 0.5) * (0.97^6)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[366:417,])), probs = 0.5) * (0.97^7)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[418:469,])), probs = 0.5) * (0.97^8)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.aids.dx <- ((quantile(unname(colSums(sim$epi$stage.time.aids.dx[2:53,])), probs = 0.5) * (0.97^0)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[54:105,])), probs = 0.5) * (0.97^1)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[106:157,])), probs = 0.5) * (0.97^2)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[158:209,])), probs = 0.5) * (0.97^3)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[210:261,])), probs = 0.5) * (0.97^4)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[262:313,])), probs = 0.5) * (0.97^5)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[314:365,])), probs = 0.5) * (0.97^6)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[366:417,])), probs = 0.5) * (0.97^7)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[418:469,])), probs = 0.5) * (0.97^8)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.aids.dx[470:521,])), probs = 0.5) * (0.97^9)))
    stage.time.aids.art <- ((quantile(unname(colSums(sim$epi$stage.time.aids.art[2:53,])), probs = 0.5) * (0.97^0)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[54:105,])), probs = 0.5) * (0.97^1)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[106:157,])), probs = 0.5) * (0.97^2)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[158:209,])), probs = 0.5) * (0.97^3)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[210:261,])), probs = 0.5) * (0.97^4)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[262:313,])), probs = 0.5) * (0.97^5)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[314:365,])), probs = 0.5) * (0.97^6)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[366:417,])), probs = 0.5) * (0.97^7)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[418:469,])), probs = 0.5) * (0.97^8)) +
                                      (quantile(unname(colSums(sim$epi$stage.time.aids.art[470:521,])), probs = 0.5) * (0.97^9)))

    # HIV/STI tests
    hivtests.nprep <- ((quantile(colSums(unname(sim$epi$hivtests.nprep[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$hivtests.nprep[470:521,])), probs = 0.5) * (0.97^9)))
    rGCsympttests <- ((quantile(colSums(unname(sim$epi$rGCsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$rGCsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    uGCsympttests <- ((quantile(colSums(unname(sim$epi$uGCsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$uGCsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    GCsympttests <- ((quantile(colSums(unname(sim$epi$GCsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$GCsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    rCTsympttests <- ((quantile(colSums(unname(sim$epi$rCTsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$rCTsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    uCTsympttests <- ((quantile(colSums(unname(sim$epi$uCTsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$uCTsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    CTsympttests <- ((quantile(colSums(unname(sim$epi$CTsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$CTsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    syphsympttests <- ((quantile(colSums(unname(sim$epi$syphsympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$syphsympttests[470:521,])), probs = 0.5) * (0.97^9)))
    rGCasympttests <- ((quantile(colSums(unname(sim$epi$rGCasympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$rGCasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    uGCasympttests <- ((quantile(colSums(unname(sim$epi$uGCasympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$uGCasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    GCasympttests <- ((quantile(colSums(unname(sim$epi$GCasympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$GCasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    rCTasympttests <- ((quantile(colSums(unname(sim$epi$rCTasympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$rCTasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    uCTasympttests <- ((quantile(colSums(unname(sim$epi$uCTasympttests[2:53, ])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$uCTasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    CTasympttests <- ((quantile(colSums(unname(sim$epi$CTasympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$CTasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    syphasympttests <- ((quantile(colSums(unname(sim$epi$syphasympttests[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$syphasympttests[470:521,])), probs = 0.5) * (0.97^9)))
    hivtests.pos <- ((quantile(colSums(unname(sim$epi$hivtests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$hivtests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    rGCasympttests.pos <- ((quantile(colSums(unname(sim$epi$rGCasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$rGCasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    uGCasympttests.pos <- ((quantile(colSums(unname(sim$epi$uGCasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$uGCasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    GCasympttests.pos <- ((quantile(colSums(unname(sim$epi$GCasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$GCasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    rCTasympttests.pos <- ((quantile(colSums(unname(sim$epi$rCTasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$rCTasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    uCTasympttests.pos <- ((quantile(colSums(unname(sim$epi$uCTasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                         (quantile(colSums(unname(sim$epi$uCTasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    CTasympttests.pos <- ((quantile(colSums(unname(sim$epi$CTasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$CTasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))
    syphasympttests.pos <- ((quantile(colSums(unname(sim$epi$syphasympttests.pos[2:53,])), probs = 0.5) * (0.97^0)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[54:105,])), probs = 0.5) * (0.97^1)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[106:157,])), probs = 0.5) * (0.97^2)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[158:209,])), probs = 0.5) * (0.97^3)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[210:261,])), probs = 0.5) * (0.97^4)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[262:313,])), probs = 0.5) * (0.97^5)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[314:365,])), probs = 0.5) * (0.97^6)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[366:417,])), probs = 0.5) * (0.97^7)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[418:469,])), probs = 0.5) * (0.97^8)) +
                          (quantile(colSums(unname(sim$epi$syphasympttests.pos[470:521,])), probs = 0.5) * (0.97^9)))

    # STI Tx
    txCT <- ((quantile(colSums(unname(sim$epi$txCT[2:53,])), probs = 0.5) * (0.97^0)) +
                              (quantile(colSums(unname(sim$epi$txCT[54:105,])), probs = 0.5) * (0.97^1)) +
                              (quantile(colSums(unname(sim$epi$txCT[106:157,])), probs = 0.5) * (0.97^2)) +
                              (quantile(colSums(unname(sim$epi$txCT[158:209,])), probs = 0.5) * (0.97^3)) +
                              (quantile(colSums(unname(sim$epi$txCT[210:261,])), probs = 0.5) * (0.97^4)) +
                              (quantile(colSums(unname(sim$epi$txCT[262:313,])), probs = 0.5) * (0.97^5)) +
                              (quantile(colSums(unname(sim$epi$txCT[314:365,])), probs = 0.5) * (0.97^6)) +
                              (quantile(colSums(unname(sim$epi$txCT[366:417,])), probs = 0.5) * (0.97^7)) +
                              (quantile(colSums(unname(sim$epi$txCT[418:469,])), probs = 0.5) * (0.97^8)) +
                              (quantile(colSums(unname(sim$epi$txCT[470:521,])), probs = 0.5) * (0.97^9)))
    txGC <- ((quantile(colSums(unname(sim$epi$txGC[2:53,])), probs = 0.5) * (0.97^0)) +
                              (quantile(colSums(unname(sim$epi$txGC[54:105,])), probs = 0.5) * (0.97^1)) +
                              (quantile(colSums(unname(sim$epi$txGC[106:157,])), probs = 0.5) * (0.97^2)) +
                              (quantile(colSums(unname(sim$epi$txGC[158:209,])), probs = 0.5) * (0.97^3)) +
                              (quantile(colSums(unname(sim$epi$txGC[210:261,])), probs = 0.5) * (0.97^4)) +
                              (quantile(colSums(unname(sim$epi$txGC[262:313,])), probs = 0.5) * (0.97^5)) +
                              (quantile(colSums(unname(sim$epi$txGC[314:365,])), probs = 0.5) * (0.97^6)) +
                              (quantile(colSums(unname(sim$epi$txGC[366:417,])), probs = 0.5) * (0.97^7)) +
                              (quantile(colSums(unname(sim$epi$txGC[418:469,])), probs = 0.5) * (0.97^8)) +
                              (quantile(colSums(unname(sim$epi$txGC[470:521,])), probs = 0.5) * (0.97^9)))
    txearlysyph <- ((quantile(colSums(unname(sim$epi$txearlysyph[2:53,])), probs = 0.5) * (0.97^0)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[54:105,])), probs = 0.5) * (0.97^1)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[106:157,])), probs = 0.5) * (0.97^2)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[158:209,])), probs = 0.5) * (0.97^3)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[210:261,])), probs = 0.5) * (0.97^4)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[262:313,])), probs = 0.5) * (0.97^5)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[314:365,])), probs = 0.5) * (0.97^6)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[366:417,])), probs = 0.5) * (0.97^7)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[418:469,])), probs = 0.5) * (0.97^8)) +
                              (quantile(colSums(unname(sim$epi$txearlysyph[470:521,])), probs = 0.5) * (0.97^9)))
    txlatesyph <- ((quantile(colSums(unname(sim$epi$txlatesyph[2:53,])), probs = 0.5) * (0.97^0)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[54:105,])), probs = 0.5) * (0.97^1)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[106:157,])), probs = 0.5) * (0.97^2)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[158:209,])), probs = 0.5) * (0.97^3)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[210:261,])), probs = 0.5) * (0.97^4)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[262:313,])), probs = 0.5) * (0.97^5)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[314:365,])), probs = 0.5) * (0.97^6)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[366:417,])), probs = 0.5) * (0.97^7)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[418:469,])), probs = 0.5) * (0.97^8)) +
                              (quantile(colSums(unname(sim$epi$txlatesyph[470:521,])), probs = 0.5) * (0.97^9)))
    txasympt <- ((quantile(colSums(unname(sim$epi$txasympt[2:53,])), probs = 0.5) * (0.97^0)) +
                              (quantile(colSums(unname(sim$epi$txasympt[54:105,])), probs = 0.5) * (0.97^1)) +
                              (quantile(colSums(unname(sim$epi$txasympt[106:157,])), probs = 0.5) * (0.97^2)) +
                              (quantile(colSums(unname(sim$epi$txasympt[158:209,])), probs = 0.5) * (0.97^3)) +
                              (quantile(colSums(unname(sim$epi$txasympt[210:261,])), probs = 0.5) * (0.97^4)) +
                              (quantile(colSums(unname(sim$epi$txasympt[262:313,])), probs = 0.5) * (0.97^5)) +
                              (quantile(colSums(unname(sim$epi$txasympt[314:365,])), probs = 0.5) * (0.97^6)) +
                              (quantile(colSums(unname(sim$epi$txasympt[366:417,])), probs = 0.5) * (0.97^7)) +
                              (quantile(colSums(unname(sim$epi$txasympt[418:469,])), probs = 0.5) * (0.97^8)) +
                              (quantile(colSums(unname(sim$epi$txasympt[470:521,])), probs = 0.5) * (0.97^9)))

    # QALY
    df$QALY[i] <- sum((time.hivneg * 1),
                      (stage.time.ar.ndx * 0.92), (stage.time.ar.dx * 0.86),
                      (stage.time.af.ndx * 0.92),  (stage.time.af.dx * 0.86),
                      (stage.time.early.chronic.ndx * 0.91), (stage.time.early.chronic.dx.yrone * 0.84),
                      (stage.time.early.chronic.dx.yrstwotolate * 0.89), (stage.time.early.chronic.art * 0.95)
                      (stage.time.late.chronic.ndx * 0.79), (stage.time.late.chronic.dx * 0.72),
                      (stage.time.late.chronic.art * 0.83),(stage.time.aids.ndx * 0.72),
                      (stage.time.aids.dx * 0.72), (stage.time.aids.art * 0.82))

    df$stiasympttests[i] <- sum((rGCasympttests), (uGCasympttests),
                             (rCTasympttests), (uCTasympttests),
                             (syphasympttests))
    df$hivtestcosts[i] <- sum((hivtests.nprep * 38.80), (hivtests.pos * 30.36))
    df$hivhealthcosts[i] <- sum((time.hivneg * 4362.67),
                                (stage.time.ar.ndx * 4667.17), (stage.time.ar.dx * 4667.17),
                                (stage.time.af.ndx * 4667.17),  (stage.time.af.dx * 4667.17),
                                (stage.time.early.chronic.ndx * 9376.9), (stage.time.early.chronic.dx.yrone * 9376.9),
                                (stage.time.early.chronic.dx.yrstwotolate * 9376.9), (stage.time.early.chronic.art * 26421.62)
                                (stage.time.late.chronic.ndx * 12597.23), (stage.time.late.chronic.dx * 12597.23),
                                (stage.time.late.chronic.art * 29641.95),(stage.time.aids.ndx * 29748.91),
                                (stage.time.aids.dx * 29748.91), (stage.time.aids.art * 33970.99))

    df$stisympttestcosts[i] <- sum((rGCsympttests * (180.64)), (uGCsympttests * (180.64)),
                                   (rCTsympttests * (180.64)), (uCTsympttests * (180.64)),
                                   (syphsympttests * (168.70)))
    df$gc.asympttestcosts[i] <- sum((GCasympttests * (180.64)))
    df$ct.asympttestcosts[i] <- sum((CTasympttests * (180.64)))
    df$syph.asympttestcosts[i] <- sum((syphasympttests * (168.70)))
    df$rect.asympttestcosts[i] <- sum((rCTasympttests * (180.64)), (rGCasympttests * (180.64)))
    df$ureth.asympttestcosts[i] <- sum((uCTasympttests * (180.64)), (uGCasympttests * (180.64)))
    df$stiasympttestcosts[i] <- sum(sum((rGCasympttests * (180.64)), (uGCasympttests * (180.64)),
                                        (rCTasympttests * (180.64)), (uCTasympttests * (180.64)),
                                        (syphasympttests * (168.70))))
    df$stitxcosts[i] <- sum((txlatesyph * 115.05), (txearlysyph * 57.52), (txGC * 51.91), (txCT * 51.91))

    df$overallcosts[i] <- sum(df$hivtestcosts[i], df$hivhealthcosts[i], df$stisympttestcosts[i], df$stiasympttestcosts[i], df$stitxcosts[i])

    df$incremover[1] <- 0

    df$incremQALYover[1] <- 0

    df$costpqaly[1] <- df$overallcosts[i] / df$QALY[i]
    df$costpoverqaly[1] <- 0


    if (i >= 2) {

    df$incremover[i] <- df$stiasympttestcosts[i] - df$stiasympttestcosts[1]
    df$incremQALYover[i] <- df$QALY[i] - df$QALY[1]
    df$costpqaly[i] <- df$overallcosts[i] / df$QALY[i]
    df$costpoverqaly[i] <- (df$incremover[i] / df$incremQALYover[i])

    }

    cat("*")

}

write.csv(df, "C:/Users/kweiss2/Documents/GitHub/stitestguidelines/analysis/STD CEA Table 1.csv")

# On laptop
#write.csv(df, "/Users/kvnweiss/stitestguidelines/analysis/Table 1.csv")
