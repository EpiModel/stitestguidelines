## STI Testing Guidelines CEA Table 2
# Varying intervals at 20% coverage

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Base STI lower-risk testing interval (364 days): n3032
## Varying STI lower-risk testing interval: 3131, 3132, 3133, 3134
## Base STI higher-risk testing interval: n3012
## Varying STI higher-risk testing interval: 3135, 3136, 3137, 3138
sims <- c(3000, 3131:3132, 3032, 3133:3134, 3135:3136, 3012, 3137:3138)

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
  time.hivneg <- ((quantile(unname(colSums(sim$epi$time.hivneg[1:52,])), probs = 0.5) * (0.97^0)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[53:104,])), probs = 0.5) * (0.97^1)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[105:156,])), probs = 0.5) * (0.97^2)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[157:208])), probs = 0.5) * (0.97^3)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[209:260,])), probs = 0.5) * (0.97^4)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[261:312,])), probs = 0.5) * (0.97^5)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[313:364,])), probs = 0.5) * (0.97^6)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[365:416,])), probs = 0.5) * (0.97^7)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[417:468,])), probs = 0.5) * (0.97^8)) +
                    (quantile(unname(colSums(sim$epi$time.hivneg[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.ar.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.ar.ndx[1:52,])), probs = 0.5) * (0.97^0)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[53:104,])), probs = 0.5) * (0.97^1)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[105:156,])), probs = 0.5) * (0.97^2)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[157:208])), probs = 0.5) * (0.97^3)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[209:260,])), probs = 0.5) * (0.97^4)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[261:312,])), probs = 0.5) * (0.97^5)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[313:364,])), probs = 0.5) * (0.97^6)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[365:416,])), probs = 0.5) * (0.97^7)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[417:468,])), probs = 0.5) * (0.97^8)) +
                          (quantile(unname(colSums(sim$epi$stage.time.ar.ndx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.ar.dx <- ((quantile(unname(colSums(sim$epi$stage.time.ar.dx[1:52,])), probs = 0.5) * (0.97^0)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[53:104,])), probs = 0.5) * (0.97^1)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[105:156,])), probs = 0.5) * (0.97^2)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[157:208])), probs = 0.5) * (0.97^3)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[209:260,])), probs = 0.5) * (0.97^4)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[261:312,])), probs = 0.5) * (0.97^5)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[313:364,])), probs = 0.5) * (0.97^6)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[365:416,])), probs = 0.5) * (0.97^7)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[417:468,])), probs = 0.5) * (0.97^8)) +
                         (quantile(unname(colSums(sim$epi$stage.time.ar.dx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.af.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.af.ndx[1:52,])), probs = 0.5) * (0.97^0)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[53:104,])), probs = 0.5) * (0.97^1)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[105:156,])), probs = 0.5) * (0.97^2)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[157:208])), probs = 0.5) * (0.97^3)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[209:260,])), probs = 0.5) * (0.97^4)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[261:312,])), probs = 0.5) * (0.97^5)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[313:364,])), probs = 0.5) * (0.97^6)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[365:416,])), probs = 0.5) * (0.97^7)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[417:468,])), probs = 0.5) * (0.97^8)) +
                          (quantile(unname(colSums(sim$epi$stage.time.af.ndx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.af.dx <- ((quantile(unname(colSums(sim$epi$stage.time.af.dx[1:52,])), probs = 0.5) * (0.97^0)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[53:104,])), probs = 0.5) * (0.97^1)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[105:156,])), probs = 0.5) * (0.97^2)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[157:208])), probs = 0.5) * (0.97^3)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[209:260,])), probs = 0.5) * (0.97^4)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[261:312,])), probs = 0.5) * (0.97^5)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[313:364,])), probs = 0.5) * (0.97^6)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[365:416,])), probs = 0.5) * (0.97^7)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[417:468,])), probs = 0.5) * (0.97^8)) +
                         (quantile(unname(colSums(sim$epi$stage.time.af.dx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.early.chronic.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[1:52,])), probs = 0.5) * (0.97^0)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[53:104,])), probs = 0.5) * (0.97^1)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[105:156,])), probs = 0.5) * (0.97^2)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[157:208])), probs = 0.5) * (0.97^3)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[209:260,])), probs = 0.5) * (0.97^4)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[261:312,])), probs = 0.5) * (0.97^5)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[313:364,])), probs = 0.5) * (0.97^6)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[365:416,])), probs = 0.5) * (0.97^7)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[417:468,])), probs = 0.5) * (0.97^8)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.ndx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.early.chronic.dx.yrone <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[1:52,])), probs = 0.5) * (0.97^0)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[53:104,])), probs = 0.5) * (0.97^1)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[105:156,])), probs = 0.5) * (0.97^2)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[157:208])), probs = 0.5) * (0.97^3)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[209:260,])), probs = 0.5) * (0.97^4)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[261:312,])), probs = 0.5) * (0.97^5)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[313:364,])), probs = 0.5) * (0.97^6)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[365:416,])), probs = 0.5) * (0.97^7)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[417:468,])), probs = 0.5) * (0.97^8)) +
                                          (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrone[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.early.chronic.dx.yrstwotolate <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[1:52,])), probs = 0.5) * (0.97^0)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[53:104,])), probs = 0.5) * (0.97^1)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[105:156,])), probs = 0.5) * (0.97^2)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[157:208])), probs = 0.5) * (0.97^3)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[209:260,])), probs = 0.5) * (0.97^4)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[261:312,])), probs = 0.5) * (0.97^5)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[313:364,])), probs = 0.5) * (0.97^6)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[365:416,])), probs = 0.5) * (0.97^7)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[417:468,])), probs = 0.5) * (0.97^8)) +
                                                 (quantile(unname(colSums(sim$epi$stage.time.early.chronic.dx.yrstwotolate[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.early.chronic.art <- ((quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[1:52,])), probs = 0.5) * (0.97^0)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[53:104,])), probs = 0.5) * (0.97^1)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[105:156,])), probs = 0.5) * (0.97^2)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[157:208])), probs = 0.5) * (0.97^3)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[209:260,])), probs = 0.5) * (0.97^4)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[261:312,])), probs = 0.5) * (0.97^5)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[313:364,])), probs = 0.5) * (0.97^6)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[365:416,])), probs = 0.5) * (0.97^7)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[417:468,])), probs = 0.5) * (0.97^8)) +
                                     (quantile(unname(colSums(sim$epi$stage.time.early.chronic.art[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.late.chronic.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[1:52,])), probs = 0.5) * (0.97^0)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[53:104,])), probs = 0.5) * (0.97^1)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[105:156,])), probs = 0.5) * (0.97^2)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[157:208])), probs = 0.5) * (0.97^3)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[209:260,])), probs = 0.5) * (0.97^4)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[261:312,])), probs = 0.5) * (0.97^5)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[313:364,])), probs = 0.5) * (0.97^6)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[365:416,])), probs = 0.5) * (0.97^7)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[417:468,])), probs = 0.5) * (0.97^8)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.ndx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.late.chronic.dx <- ((quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[1:52,])), probs = 0.5) * (0.97^0)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[53:104,])), probs = 0.5) * (0.97^1)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[105:156,])), probs = 0.5) * (0.97^2)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[157:208])), probs = 0.5) * (0.97^3)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[209:260,])), probs = 0.5) * (0.97^4)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[261:312,])), probs = 0.5) * (0.97^5)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[313:364,])), probs = 0.5) * (0.97^6)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[365:416,])), probs = 0.5) * (0.97^7)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[417:468,])), probs = 0.5) * (0.97^8)) +
                                   (quantile(unname(colSums(sim$epi$stage.time.late.chronic.dx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.late.chronic.art <- ((quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[1:52,])), probs = 0.5) * (0.97^0)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[53:104,])), probs = 0.5) * (0.97^1)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[105:156,])), probs = 0.5) * (0.97^2)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[157:208])), probs = 0.5) * (0.97^3)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[209:260,])), probs = 0.5) * (0.97^4)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[261:312,])), probs = 0.5) * (0.97^5)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[313:364,])), probs = 0.5) * (0.97^6)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[365:416,])), probs = 0.5) * (0.97^7)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[417:468,])), probs = 0.5) * (0.97^8)) +
                                    (quantile(unname(colSums(sim$epi$stage.time.late.chronic.art[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.aids.ndx <- ((quantile(unname(colSums(sim$epi$stage.time.aids.ndx[1:52,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[53:104,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[105:156,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[157:208])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[209:260,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[261:312,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[313:364,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[365:416,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[417:468,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.ndx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.aids.dx <- ((quantile(unname(colSums(sim$epi$stage.time.aids.dx[1:52,])), probs = 0.5) * (0.97^0)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[53:104,])), probs = 0.5) * (0.97^1)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[105:156,])), probs = 0.5) * (0.97^2)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[157:208])), probs = 0.5) * (0.97^3)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[209:260,])), probs = 0.5) * (0.97^4)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[261:312,])), probs = 0.5) * (0.97^5)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[313:364,])), probs = 0.5) * (0.97^6)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[365:416,])), probs = 0.5) * (0.97^7)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[417:468,])), probs = 0.5) * (0.97^8)) +
                           (quantile(unname(colSums(sim$epi$stage.time.aids.dx[469:520,])), probs = 0.5) * (0.97^9)))
  stage.time.aids.art <- ((quantile(unname(colSums(sim$epi$stage.time.aids.art[1:52,])), probs = 0.5) * (0.97^0)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[53:104,])), probs = 0.5) * (0.97^1)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[105:156,])), probs = 0.5) * (0.97^2)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[157:208])), probs = 0.5) * (0.97^3)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[209:260,])), probs = 0.5) * (0.97^4)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[261:312,])), probs = 0.5) * (0.97^5)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[313:364,])), probs = 0.5) * (0.97^6)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[365:416,])), probs = 0.5) * (0.97^7)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[417:468,])), probs = 0.5) * (0.97^8)) +
                            (quantile(unname(colSums(sim$epi$stage.time.aids.art[469:520,])), probs = 0.5) * (0.97^9)))

  # HIV/STI tests
  hivtests.nprep <- ((quantile(colSums(unname(sim$epi$hivtests.nprep[1:52, ])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[53:104, ])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[105:156, ])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[157:208, ])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[209:260, ])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[261:312, ])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[313:364, ])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[365:416, ])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[417:468, ])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$hivtests.nprep[469:520, ])), probs = 0.5) * (0.97^9)))
  rGCsympttests <- ((quantile(colSums(unname(sim$epi$rGCsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                      (quantile(colSums(unname(sim$epi$rGCsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  uGCsympttests <- ((quantile(colSums(unname(sim$epi$uGCsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                      (quantile(colSums(unname(sim$epi$uGCsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  GCsympttests <- ((quantile(colSums(unname(sim$epi$GCsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                     (quantile(colSums(unname(sim$epi$GCsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  rCTsympttests <- ((quantile(colSums(unname(sim$epi$rCTsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                      (quantile(colSums(unname(sim$epi$rCTsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  uCTsympttests <- ((quantile(colSums(unname(sim$epi$uCTsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                      (quantile(colSums(unname(sim$epi$uCTsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  CTsympttests <- ((quantile(colSums(unname(sim$epi$CTsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                     (quantile(colSums(unname(sim$epi$CTsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  syphsympttests <- ((quantile(colSums(unname(sim$epi$syphsympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$syphsympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  rGCasympttests <- ((quantile(colSums(unname(sim$epi$rGCasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$rGCasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  uGCasympttests <- ((quantile(colSums(unname(sim$epi$uGCasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$uGCasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  GCasympttests <- ((quantile(colSums(unname(sim$epi$GCasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                      (quantile(colSums(unname(sim$epi$GCasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  rCTasympttests <- ((quantile(colSums(unname(sim$epi$rCTasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$rCTasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  uCTasympttests <- ((quantile(colSums(unname(sim$epi$uCTasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                       (quantile(colSums(unname(sim$epi$uCTasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  CTasympttests <- ((quantile(colSums(unname(sim$epi$CTasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                      (quantile(colSums(unname(sim$epi$CTasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  syphasympttests <- ((quantile(colSums(unname(sim$epi$syphasympttests[1:52, ])), probs = 0.5) * (0.97^0)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[53:104, ])), probs = 0.5) * (0.97^1)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[105:156, ])), probs = 0.5) * (0.97^2)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[157:208, ])), probs = 0.5) * (0.97^3)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[209:260, ])), probs = 0.5) * (0.97^4)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[261:312, ])), probs = 0.5) * (0.97^5)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[313:364, ])), probs = 0.5) * (0.97^6)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[365:416, ])), probs = 0.5) * (0.97^7)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[417:468, ])), probs = 0.5) * (0.97^8)) +
                        (quantile(colSums(unname(sim$epi$syphasympttests[469:520, ])), probs = 0.5) * (0.97^9)))
  hivtests.pos <- ((quantile(colSums(unname(sim$epi$hivtests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                     (quantile(colSums(unname(sim$epi$hivtests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  rGCasympttests.pos <- ((quantile(colSums(unname(sim$epi$rGCasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                               (quantile(colSums(unname(sim$epi$rGCasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  uGCasympttests.pos <- ((quantile(colSums(unname(sim$epi$uGCasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                           (quantile(colSums(unname(sim$epi$uGCasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  GCasympttests.pos <- ((quantile(colSums(unname(sim$epi$GCasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                          (quantile(colSums(unname(sim$epi$GCasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  rCTasympttests.pos <- ((quantile(colSums(unname(sim$epi$rCTasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                           (quantile(colSums(unname(sim$epi$rCTasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  uCTasympttests.pos <- ((quantile(colSums(unname(sim$epi$uCTasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                           (quantile(colSums(unname(sim$epi$uCTasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  CTasympttests.pos <- ((quantile(colSums(unname(sim$epi$CTasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                          (quantile(colSums(unname(sim$epi$CTasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))
  syphasympttests.pos <- ((quantile(colSums(unname(sim$epi$syphasympttests.pos[1:52, ])), probs = 0.5) * (0.97^0)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[53:104, ])), probs = 0.5) * (0.97^1)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[105:156, ])), probs = 0.5) * (0.97^2)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[157:208, ])), probs = 0.5) * (0.97^3)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[209:260, ])), probs = 0.5) * (0.97^4)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[261:312, ])), probs = 0.5) * (0.97^5)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[313:364, ])), probs = 0.5) * (0.97^6)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[365:416, ])), probs = 0.5) * (0.97^7)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[417:468, ])), probs = 0.5) * (0.97^8)) +
                            (quantile(colSums(unname(sim$epi$syphasympttests.pos[469:520, ])), probs = 0.5) * (0.97^9)))

  # STI Tx
  txCT <- ((quantile(colSums(unname(sim$epi$txCT[1:52, ])), probs = 0.5) * (0.97^0)) +
             (quantile(colSums(unname(sim$epi$txCT[53:104, ])), probs = 0.5) * (0.97^1)) +
             (quantile(colSums(unname(sim$epi$txCT[105:156, ])), probs = 0.5) * (0.97^2)) +
             (quantile(colSums(unname(sim$epi$txCT[157:208, ])), probs = 0.5) * (0.97^3)) +
             (quantile(colSums(unname(sim$epi$txCT[209:260, ])), probs = 0.5) * (0.97^4)) +
             (quantile(colSums(unname(sim$epi$txCT[261:312, ])), probs = 0.5) * (0.97^5)) +
             (quantile(colSums(unname(sim$epi$txCT[313:364, ])), probs = 0.5) * (0.97^6)) +
             (quantile(colSums(unname(sim$epi$txCT[365:416, ])), probs = 0.5) * (0.97^7)) +
             (quantile(colSums(unname(sim$epi$txCT[417:468, ])), probs = 0.5) * (0.97^8)) +
             (quantile(colSums(unname(sim$epi$txCT[469:520, ])), probs = 0.5) * (0.97^9)))
  txGC <- ((quantile(colSums(unname(sim$epi$txGC[1:52, ])), probs = 0.5) * (0.97^0)) +
             (quantile(colSums(unname(sim$epi$txGC[53:104, ])), probs = 0.5) * (0.97^1)) +
             (quantile(colSums(unname(sim$epi$txGC[105:156, ])), probs = 0.5) * (0.97^2)) +
             (quantile(colSums(unname(sim$epi$txGC[157:208, ])), probs = 0.5) * (0.97^3)) +
             (quantile(colSums(unname(sim$epi$txGC[209:260, ])), probs = 0.5) * (0.97^4)) +
             (quantile(colSums(unname(sim$epi$txGC[261:312, ])), probs = 0.5) * (0.97^5)) +
             (quantile(colSums(unname(sim$epi$txGC[313:364, ])), probs = 0.5) * (0.97^6)) +
             (quantile(colSums(unname(sim$epi$txGC[365:416, ])), probs = 0.5) * (0.97^7)) +
             (quantile(colSums(unname(sim$epi$txGC[417:468, ])), probs = 0.5) * (0.97^8)) +
             (quantile(colSums(unname(sim$epi$txGC[469:520, ])), probs = 0.5) * (0.97^9)))
  txearlysyph <- ((quantile(colSums(unname(sim$epi$txearlysyph[1:52, ])), probs = 0.5) * (0.97^0)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[53:104, ])), probs = 0.5) * (0.97^1)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[105:156, ])), probs = 0.5) * (0.97^2)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[157:208, ])), probs = 0.5) * (0.97^3)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[209:260, ])), probs = 0.5) * (0.97^4)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[261:312, ])), probs = 0.5) * (0.97^5)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[313:364, ])), probs = 0.5) * (0.97^6)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[365:416, ])), probs = 0.5) * (0.97^7)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[417:468, ])), probs = 0.5) * (0.97^8)) +
                    (quantile(colSums(unname(sim$epi$txearlysyph[469:520, ])), probs = 0.5) * (0.97^9)))
  txlatesyph <- ((quantile(colSums(unname(sim$epi$txlatesyph[1:52, ])), probs = 0.5) * (0.97^0)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[53:104, ])), probs = 0.5) * (0.97^1)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[105:156, ])), probs = 0.5) * (0.97^2)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[157:208, ])), probs = 0.5) * (0.97^3)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[209:260, ])), probs = 0.5) * (0.97^4)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[261:312, ])), probs = 0.5) * (0.97^5)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[313:364, ])), probs = 0.5) * (0.97^6)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[365:416, ])), probs = 0.5) * (0.97^7)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[417:468, ])), probs = 0.5) * (0.97^8)) +
                   (quantile(colSums(unname(sim$epi$txlatesyph[469:520, ])), probs = 0.5) * (0.97^9)))
  txasympt <- ((quantile(colSums(unname(sim$epi$txasympt[1:52, ])), probs = 0.5) * (0.97^0)) +
                 (quantile(colSums(unname(sim$epi$txasympt[53:104, ])), probs = 0.5) * (0.97^1)) +
                 (quantile(colSums(unname(sim$epi$txasympt[105:156, ])), probs = 0.5) * (0.97^2)) +
                 (quantile(colSums(unname(sim$epi$txasympt[157:208, ])), probs = 0.5) * (0.97^3)) +
                 (quantile(colSums(unname(sim$epi$txasympt[209:260, ])), probs = 0.5) * (0.97^4)) +
                 (quantile(colSums(unname(sim$epi$txasympt[261:312, ])), probs = 0.5) * (0.97^5)) +
                 (quantile(colSums(unname(sim$epi$txasympt[313:364, ])), probs = 0.5) * (0.97^6)) +
                 (quantile(colSums(unname(sim$epi$txasympt[365:416, ])), probs = 0.5) * (0.97^7)) +
                 (quantile(colSums(unname(sim$epi$txasympt[417:468, ])), probs = 0.5) * (0.97^8)) +
                 (quantile(colSums(unname(sim$epi$txasympt[469:520, ])), probs = 0.5) * (0.97^9)))

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

write.csv(df, "C:/Users/kweiss2/Documents/GitHub/stitestguidelines/analysis/STD CEA Table 2.csv")
