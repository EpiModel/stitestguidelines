library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")


sim <- truncate_sim(sim, at = 2600)

# Take value at end of simulation
time.hivneg <- as.numeric(sim$epi$time.hivneg[520, ]) / 52
round(quantile(time.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

time.on.prep <- as.numeric(sim$epi$time.on.prep[520, ]) / 52
round(quantile(time.on.prep, probs = c(0.025, 0.5, 0.975)), 3)
 
time.off.prep <- as.numeric(sim$epi$time.off.prep[520, ]) / 52
round(quantile(time.off.prep, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar <- as.numeric(sim$epi$stage.time.ar[520, ]) / 52
round(quantile(stage.time.ar, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af <- as.numeric(sim$epi$stage.time.af[520, ]) / 52
round(quantile(stage.time.af, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic <- as.numeric(sim$epi$stage.time.chronic[520, ]) / 52
round(quantile(stage.time.chronic, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids <- as.numeric(sim$epi$stage.time.aids[520,]) / 52
round(quantile(stage.time.aids, probs = c(0.025, 0.5, 0.975)), 3)

totalhivtests <- as.numeric(sim$epi$totalhivtests[520, ])
round(quantile(totalhivtests, probs = c(0.025, 0.5, 0.975)), 3)

totalhivtests.prep <- as.numeric(sim$epi$totalhivtests.prep[520, ])
round(quantile(totalhivtests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totalrGCsympttests <- as.numeric(sim$epi$totalrGCsympttests[520, ])
round(quantile(totalrGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCsympttests <- as.numeric(sim$epi$totaluGCsympttests[520, ])
round(quantile(totaluGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCsympttests <- as.numeric(sim$epi$totalGCsympttests[520, ])
round(quantile(totalGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTsympttests <- as.numeric(sim$epi$totalrCTsympttests[520, ])
round(quantile(totalrCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTsympttests <- as.numeric(sim$epi$totaluCTsympttests[520, ])
round(quantile(totaluCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTsympttests <- as.numeric(sim$epi$totalCTsympttests[520, ])
round(quantile(totalCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphsympttests <- as.numeric(sim$epi$totalsyphsympttests[520, ])
round(quantile(totalsyphsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrGCasympttests <- as.numeric(sim$epi$totalrGCasympttests[520, ])
round(quantile(totalrGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCasympttests <- as.numeric(sim$epi$totaluGCasympttests[520, ])
round(quantile(totaluGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCasympttests <- as.numeric(sim$epi$totalGCasympttests[520, ])
round(quantile(totalGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTasympttests <- as.numeric(sim$epi$totalrCTasympttests[520, ])
round(quantile(totalrCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTasympttests <- as.numeric(sim$epi$totaluCTasympttests[520, ])
round(quantile(totaluCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTasympttests <- as.numeric(sim$epi$totalCTasympttests[520, ])
round(quantile(totalCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphasympttests <- as.numeric(sim$epi$totalsyphasympttests[520, ])
round(quantile(totalsyphasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrGCasympttests.prep <- as.numeric(sim$epi$totalrGCasympttests.prep[520, ])
round(quantile(totalrGCasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCasympttests.prep <- as.numeric(sim$epi$totaluGCasympttests.prep[520, ])
round(quantile(totaluGCasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totalGCasympttests.prep <- as.numeric(sim$epi$totalGCasympttests.prep[520, ])
round(quantile(totalGCasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTasympttests.prep <- as.numeric(sim$epi$totalrCTasympttests.prep[520, ])
round(quantile(totalrCTasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTasympttests.prep <- as.numeric(sim$epi$totaluCTasympttests.prep[520, ])
round(quantile(totaluCTasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totalCTasympttests.prep <- as.numeric(sim$epi$totalCTasympttests.prep[520, ])
round(quantile(totalCTasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphasympttests.prep <- as.numeric(sim$epi$totalsyphasympttests.prep[520, ])
round(quantile(totalsyphasympttests.prep, probs = c(0.025, 0.5, 0.975)), 3)

# # Summary of 500 FU sims
# par(mfrow = c(1, 1), oma = c(0,0,2,0))
# plot(sim, y = "time.off.prep", ylab = "Time Spent in State", xlab = "Simulation Time", mean.col = "blue", qnts.col = "blue")
# plot(sim, y = "time.on.prep", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
# title("PrEP Time", outer = TRUE)
# legend("topleft", c("Off PrEP", "On PrEP"), col = c("blue", "red"), lty = c(1, 1))
# 
# par(mfrow = c(1, 2), oma = c(0,0,2,0))
# plot(sim, y = "time.hivneg", ylab = "Time Spent in State", mean.col = "blue", qnts.col = "blue")
# title("Time spent HIV-negative")
# plot(sim, y = "stage.time.chronic", ylab = "Time Spent in State", mean.col = "gray", qnts.col = "gray")
# plot(sim, y = "stage.time.ar", ylab = "Time Spent in State", add = TRUE, mean.col = "blue", qnts.col = "blue")
# plot(sim, y = "stage.time.af", ylab = "Time Spent in State", add = TRUE, mean.col = "red", qnts.col = "red")
# plot(sim, y = "stage.time.aids", ylab = "Time Spent in State", add = TRUE, mean.col = "purple", qnts.col = "purple")
# title("Time spent per HIV stage")
# legend("topleft", c("Acute Rising", "Acute Falling", "Chronic", "AIDS"), col = c("blue", "red", "gray", "purple"), lty = c(1, 1))
# 
# par(mfrow = c(1, 2), oma = c(0,0,2,0))
# plot(sim, y = "hivtests", ylab = "Number", mean.col = "blue", qnts.col = "blue")
# title("HIV Tests per time step")
# plot(sim, y = "totalhivtests", ylab = "Number", mean.col = "blue", qnts.col = "blue")
# title("Total HIV Tests")

df <- as.data.frame(cbind(rbind(quantile(time.hivneg, probs = 0.25), quantile(time.on.prep, probs = 0.25), 
                 # quantile(time.off.prep, probs = 0.25), 
                  quantile(stage.time.ar, probs = 0.25),
                  quantile(stage.time.af, probs = 0.25), quantile(stage.time.chronic, probs = 0.25),
                  quantile(stage.time.aids, probs = 0.25)),
            rbind(quantile(time.hivneg, probs = 0.5), quantile(time.on.prep, probs = 0.5), 
                  #quantile(time.off.prep, probs = 0.5),
                  quantile(stage.time.ar, probs = 0.5),
                  quantile(stage.time.af, probs = 0.5), quantile(stage.time.chronic, probs = 0.5), 
                  quantile(stage.time.aids, probs = 0.5)),
            rbind(quantile(time.hivneg, probs = 0.75), quantile(time.on.prep, probs = 0.75), 
                  #quantile(time.off.prep, probs = 0.75), 
                  quantile(stage.time.ar, probs = 0.75),
                  quantile(stage.time.af, probs = 0.75), quantile(stage.time.chronic, probs = 0.75), 
                  quantile(stage.time.aids, probs = 0.75))))
rownames(df) <- c("time.hivneg", "time.on.prep", #"time.off.prep", 
                  "stage.time.ar", "stage.time.af", "stage.time.chronic", "stage.time.aids")
colnames(df) <- c("Twentyfive", "Fifty", "Seventyfive")

df$weight <- c(1, 1, 0.75, 0.75, 0.5, 0.1)
df$Twentyfivevalue <- df$weight*df$Twentyfive
df$Fiftyvalue <- df$weight*df$Fifty
df$Seventyfivevalue <- df$weight*df$Seventyfive


# For first
QALY <- cbind(colSums(df[5]) / 52, colSums(df[6]), colSums(df[7]),
              quantile(totalhivtests, probs = 0.25), quantile(totalhivtests, probs = 0.5), quantile(totalhivtests, probs = 0.75),
              quantile(totalhivtests.prep, probs = 0.25), quantile(totalhivtests.prep, probs = 0.5), quantile(totalhivtests.prep, probs = 0.75),
              quantile(time.on.prep, probs = 0.25), quantile(time.on.prep, probs = 0.5), quantile(time.on.prep, probs = 0.75),
              quantile(totalrGCasympttests.prep, probs = 0.25), quantile(totalrGCasympttests.prep, probs = 0.5), quantile(totalrGCasympttests.prep, probs = 0.75),
              quantile(totaluGCasympttests.prep, probs = 0.25), quantile(totaluGCasympttests.prep, probs = 0.5), quantile(totaluGCasympttests.prep, probs = 0.75),
              quantile(totalGCasympttests.prep, probs = 0.25), quantile(totalGCasympttests.prep, probs = 0.5), quantile(totalGCasympttests.prep, probs = 0.75),
              quantile(totalrCTasympttests.prep, probs = 0.25), quantile(totalrGCasympttests.prep, probs = 0.5), quantile(totalrGCasympttests.prep, probs = 0.75),
              quantile(totaluCTasympttests.prep, probs = 0.25), quantile(totaluCTasympttests.prep, probs = 0.5), quantile(totaluCTasympttests.prep, probs = 0.75),
              quantile(totalCTasympttests.prep, probs = 0.25), quantile(totalCTasympttests.prep, probs = 0.5), quantile(totalCTasympttests.prep, probs = 0.75),
              quantile(totalsyphasympttests.prep, probs = 0.25), quantile(totalsyphasympttests.prep, probs = 0.5), quantile(totalsyphasympttests.prep, probs = 0.75))

QALY2 <- QALY
QALY2

# For others
sim <- truncate_sim(sim, at = 2600)

# Take value at end of simulation
time.hivneg <- as.numeric(sim$epi$time.hivneg[520, ])
round(quantile(time.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

time.on.prep <- as.numeric(sim$epi$time.on.prep[520, ])
round(quantile(time.on.prep, probs = c(0.025, 0.5, 0.975)), 3)

time.off.prep <- as.numeric(sim$epi$time.off.prep[520, ])
round(quantile(time.off.prep, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar <- as.numeric(sim$epi$stage.time.ar[520, ])
round(quantile(stage.time.ar, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af <- as.numeric(sim$epi$stage.time.af[520, ])
round(quantile(stage.time.af, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic <- as.numeric(sim$epi$stage.time.chronic[520, ])
round(quantile(stage.time.chronic, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids <- as.numeric(sim$epi$stage.time.aids[520,])
round(quantile(stage.time.aids, probs = c(0.025, 0.5, 0.975)), 3)

totalhivtests <- as.numeric(sim$epi$totalhivtests[520, ])
round(quantile(totalhivtests, probs = c(0.025, 0.5, 0.975)), 3)


# # Summary of 500 FU sims
# par(mfrow = c(1, 1), oma = c(0,0,2,0))
# plot(sim, y = "time.off.prep", ylab = "Time Spent in State", xlab = "Simulation Time", mean.col = "blue", qnts.col = "blue")
# plot(sim, y = "time.on.prep", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
# title("PrEP Time", outer = TRUE)
# legend("topleft", c("Off PrEP", "On PrEP"), col = c("blue", "red"), lty = c(1, 1))
# 
# par(mfrow = c(1, 2), oma = c(0,0,2,0))
# plot(sim, y = "time.hivneg", ylab = "Time Spent in State", mean.col = "blue", qnts.col = "blue")
# title("Time spent HIV-negative")
# plot(sim, y = "stage.time.chronic", ylab = "Time Spent in State", mean.col = "gray", qnts.col = "gray")
# plot(sim, y = "stage.time.ar", ylab = "Time Spent in State", add = TRUE, mean.col = "blue", qnts.col = "blue")
# plot(sim, y = "stage.time.af", ylab = "Time Spent in State", add = TRUE, mean.col = "red", qnts.col = "red")
# plot(sim, y = "stage.time.aids", ylab = "Time Spent in State", add = TRUE, mean.col = "purple", qnts.col = "purple")
# title("Time spent per HIV stage")
# legend("topleft", c("Acute Rising", "Acute Falling", "Chronic", "AIDS"), col = c("blue", "red", "gray", "purple"), lty = c(1, 1))
# 
# par(mfrow = c(1, 2), oma = c(0,0,2,0))
# plot(sim, y = "hivtests", ylab = "Number", mean.col = "blue", qnts.col = "blue")
# title("HIV Tests per time step")
# plot(sim, y = "totalhivtests", ylab = "Number", mean.col = "blue", qnts.col = "blue")
# title("Total HIV Tests")

df <- as.data.frame(cbind(rbind(quantile(time.hivneg, probs = 0.25), quantile(time.on.prep, probs = 0.25), 
                                # quantile(time.off.prep, probs = 0.25), 
                                quantile(stage.time.ar, probs = 0.25),
                                quantile(stage.time.af, probs = 0.25), quantile(stage.time.chronic, probs = 0.25),
                                quantile(stage.time.aids, probs = 0.25)),
                          rbind(quantile(time.hivneg, probs = 0.5), quantile(time.on.prep, probs = 0.5), 
                                #quantile(time.off.prep, probs = 0.5),
                                quantile(stage.time.ar, probs = 0.5),
                                quantile(stage.time.af, probs = 0.5), quantile(stage.time.chronic, probs = 0.5), 
                                quantile(stage.time.aids, probs = 0.5)),
                          rbind(quantile(time.hivneg, probs = 0.75), quantile(time.on.prep, probs = 0.75), 
                                #quantile(time.off.prep, probs = 0.75), 
                                quantile(stage.time.ar, probs = 0.75),
                                quantile(stage.time.af, probs = 0.75), quantile(stage.time.chronic, probs = 0.75), 
                                quantile(stage.time.aids, probs = 0.75))))
rownames(df) <- c("time.hivneg", "time.on.prep", #"time.off.prep", 
                  "stage.time.ar", "stage.time.af", "stage.time.chronic", "stage.time.aids")
colnames(df) <- c("Twentyfive", "Fifty", "Seventyfive")

df$weight <- c(1, 1, 0.75, 0.75, 0.5, 0.1)
df$Twentyfivevalue <- df$weight*df$Twentyfive
df$Fiftyvalue <- df$weight*df$Fifty
df$Seventyfivevalue <- df$weight*df$Seventyfive

QALY <- cbind(colSums(df[5]) / 52, colSums(df[6]) / 52, colSums(df[7]) / 52,
              quantile(totalhivtests, probs = 0.25), quantile(totalhivtests, probs = 0.5), quantile(totalhivtests, probs = 0.75), 
              quantile(time.on.prep, probs = 0.25), quantile(time.on.prep, probs = 0.5), quantile(time.on.prep, probs = 0.75))
QALY2 <- rbind(QALY2, QALY)


# Create Table for Output
rownames(QALY2) <- c("Cov = 0%", "Cov = 10%", "Cov = 20%", "Cov = 30%", "Cov = 40%", "Cov = 50%",
                     "Cov = 60%", "Cov = 70%", "Cov = 80%", "Cov = 90%", "Cov = 100%")
colnames(QALY2) <- c("25% QALY", "50% QALY", "75% QALY", "25% HIV Tests", "50% HIV Tests", "75% HIV Tests", "25% Years on PrEP", "50% Years on PrEP", "75% Years on PrEP")
QALY2
