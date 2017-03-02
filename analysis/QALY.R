library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# fn <- list.files("data/", pattern = "n[3-4][0-9][0-9][0-9].rda")
# for (i in fn) {
#     load(i)
#     sim <- truncate_sim(sim, at = 2600)
#     save(sim, file = i, compress = FALSE)
#     cat("*")
# }
# 
# 
# sims <- 3000:3011
# cov <- rep(NA, length(sims))
# rc <- rep(NA, length(sims))
# hr.gc <- rep(NA, length(sims))
# hr.ct <- rep(NA, length(sims))
# pia.gc <- rep(NA, length(sims))
# nnt.gc <- rep(NA, length(sims))
# pia.ct <- rep(NA, length(sims))
# nnt.ct <- rep(NA, length(sims))
# time.hivneg <- rep(NA, length(sims))
# stage.time.ar <- rep(NA, length(sims))
# stage.time.af <- rep(NA, length(sims))
# stage.time.chronic <- rep(NA, length(sims))
# stage.time.aids <- rep(NA, length(sims))
# 
# df <- data.frame(sims, cov, rc, hr.gc, hr.ct, pia.gc, pia.ct, nnt.gc, nnt.ct)
# 
# load("data/sim.n100.rda")
# sim.base <- sim
# haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
# ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
# incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
# 
# haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
# ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
# incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
# 
# for (i in seq_along(sims)) {
#     fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
#     load(fn)
#     df$cov[i] <- sim$param$prep.coverage
#     df$rc[i] <- sim$param$rcomp.prob
#     
#     # HR
#     num.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
#     denom.gc <- unname(colMeans(tail(sim.base$epi$ir100.gc, 52)))
#     vec.hr.gc <- num.gc/denom.gc
#     vec.hr.gc <- vec.hr.gc[vec.hr.gc < Inf]
#     df$hr.gc[i] <- median(vec.hr.gc, na.rm = TRUE)
#     
#     num.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
#     denom.ct <- unname(colMeans(tail(sim.base$epi$ir100.ct, 52)))
#     vec.hr.ct <- num.ct/denom.ct
#     vec.hr.ct <- vec.hr.ct[vec.hr.ct < Inf]
#     df$hr.ct[i] <- median(vec.hr.ct, na.rm = TRUE)
#     
#     # PIA
#     ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
#     vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
#     vec.pia.gc <- vec.nia.gc/ir.base.gc
#     vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
#     df$pia.gc[i] <- median(vec.pia.gc, na.rm = TRUE)
#     
#     ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
#     vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
#     vec.pia.ct <- vec.nia.ct/ir.base.ct
#     vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
#     df$pia.ct[i] <- median(vec.pia.ct, na.rm = TRUE)
#     
#     # NNT
#     py.on.prep <- unname(colSums(sim$epi$prepCurr))/52
#     vec.nnt.gc <- py.on.prep / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
#     vec.nnt.ct <- py.on.prep / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))
#     
#     df$nnt.gc[i] <- median(vec.nnt.gc, na.rm = TRUE)
#     df$nnt.ct[i] <- median(vec.nnt.ct, na.rm = TRUE)
#     cat("*")
# }
# df

# sim <- merge_simfiles(simno = 3000, indir = "data/", ftype = "max")

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

totalhivtests.prep <- as.numeric(sim$epi$totalhivtests.prep[520, ])
round(quantile(totalhivtests.prep, probs = c(0.025, 0.5, 0.975)), 3)

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
QALY <- cbind(colSums(df[5]) / 52, colSums(df[6]) / 52, colSums(df[7]) / 52,
              quantile(totalhivtests, probs = 0.25), quantile(totalhivtests, probs = 0.5), quantile(totalhivtests, probs = 0.75),
              quantile(totalhivtests.prep, probs = 0.25), quantile(totalhivtests.prep, probs = 0.5), quantile(totalhivtests.prep, probs = 0.75),
              quantile(time.on.prep, probs = 0.25) / 52 , quantile(time.on.prep, probs = 0.5) / 52, quantile(time.on.prep, probs = 0.75) / 52)
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
