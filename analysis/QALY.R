library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

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


# Summary of 500 FU sims
par(mfrow = c(1, 1), oma = c(0,0,2,0))
plot(sim, y = "time.off.prep", ylab = "Time Spent in State", xlab = "Simulation Time", mean.col = "blue", qnts.col = "blue")
plot(sim, y = "time.on.prep", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
title("PrEP Time", outer = TRUE)
legend("topleft", c("Off PrEP", "On PrEP"), col = c("blue", "red"), lty = c(1, 1))

par(mfrow = c(1, 2), oma = c(0,0,2,0))
plot(sim, y = "time.hivneg", ylab = "Time Spent in State", mean.col = "blue", qnts.col = "blue")
title("Time spent HIV-negative")
plot(sim, y = "stage.time.chronic", ylab = "Time Spent in State", mean.col = "gray", qnts.col = "gray")
plot(sim, y = "stage.time.ar", ylab = "Time Spent in State", add = TRUE, mean.col = "blue", qnts.col = "blue")
plot(sim, y = "stage.time.af", ylab = "Time Spent in State", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "stage.time.aids", ylab = "Time Spent in State", add = TRUE, mean.col = "purple", qnts.col = "purple")
title("Time spent per HIV stage")
legend("topleft", c("Acute Rising", "Acute Falling", "Chronic", "AIDS"), col = c("blue", "red", "gray", "purple"), lty = c(1, 1))

par(mfrow = c(1, 2), oma = c(0,0,2,0))
plot(sim, y = "hivtests", ylab = "Number", mean.col = "blue", qnts.col = "blue")
title("HIV Tests per time step")
plot(sim, y = "totalhivtests", ylab = "Number", mean.col = "blue", qnts.col = "blue")
title("Total HIV Tests")

df <- as.data.frame(cbind(rbind(quantile(time.hivneg, probs = 0.25), quantile(time.on.prep, probs = 0.25), 
                  quantile(time.off.prep, probs = 0.25), quantile(stage.time.ar, probs = 0.25),
                  quantile(stage.time.af, probs = 0.25), quantile(stage.time.chronic, probs = 0.25),
                  quantile(stage.time.aids, probs = 0.25)),
            rbind(quantile(time.hivneg, probs = 0.5), quantile(time.on.prep, probs = 0.5), 
                  quantile(time.off.prep, probs = 0.5), quantile(stage.time.ar, probs = 0.5),
                  quantile(stage.time.af, probs = 0.5), quantile(stage.time.chronic, probs = 0.5), 
                  quantile(stage.time.aids, probs = 0.5)),
            rbind(quantile(time.hivneg, probs = 0.75), quantile(time.on.prep, probs = 0.75), 
                  quantile(time.off.prep, probs = 0.75), quantile(stage.time.ar, probs = 0.75),
                  quantile(stage.time.af, probs = 0.75), quantile(stage.time.chronic, probs = 0.75), 
                  quantile(stage.time.aids, probs = 0.75))))
rownames(df) <- c("time.hivneg", "time.on.prep", "time.off.prep", "stage.time.ar", "stage.time.af", "stage.time.chronic", "stage.time.aids")
colnames(df) <- c("Twentyfive", "Fifty", "Seventyfive")
df$weight <- c(1, 1, 1, 0.75, 0.75, 0.5, 0.1)
df$Twentyfivevalue <- df$weight*df$Twentyfive
df$Fiftyvalue <- df$weight*df$Fifty
df$Seventyfivevalue <- df$weight*df$Seventyfive
QALY <- colSums(df[5:7])
QALY
