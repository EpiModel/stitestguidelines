## Incidence

ir100.gc <- as.numeric(sim$epi$ir100.gc[2600, ])
round(quantile(ir100.gc, probs = c(0.025, 0.5, 0.975)), 3)

ir100.ct <- as.numeric(sim$epi$ir100.ct[2600, ])
round(quantile(ir100.ct, probs = c(0.025, 0.5, 0.975)), 3)

ir100 <- as.numeric(sim$epi$ir100[2600, ])
round(quantile(ir100, probs = c(0.025, 0.5, 0.975)), 3)

ir100.syph <- as.numeric(sim$epi$ir100.syph[2600, ])
round(quantile(ir100.syph, probs = c(0.025, 0.5, 0.975)), 3)

i.prev <- as.numeric(sim$epi$i.prev[2600, ])
round(quantile(i.prev, probs = c(0.025, 0.5, 0.975)), 3)

# Syphilis stage prevalence and HIV-syphilis overlap

prev.primsecosyph.hivpos <- as.numeric(sim$epi$prev.primsecosyph.hivpos[2600, ])
round(quantile(prev.primsecosyph.hivpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.primsecosyph.hivneg <- as.numeric(sim$epi$prev.primsecosyph.hivneg[2600, ])
round(quantile(prev.primsecosyph.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph.hivpos <- as.numeric(sim$epi$prev.syph.hivpos[2600, ])
round(quantile(prev.syph.hivpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph.hivneg <- as.numeric(sim$epi$prev.syph.hivneg[2600, ])
round(quantile(prev.syph.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph <- as.numeric(sim$epi$prev.syph[2600, ])
round(quantile(prev.syph, probs = c(0.025, 0.5, 0.975)), 3)

prev.primsecosyph <- as.numeric(sim$epi$prev.primsecosyph[2600, ])
round(quantile(prev.primsecosyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.hiv.syphpos <- as.numeric(sim$epi$prev.hiv.syphpos[2600, ])
round(quantile(prev.hiv.syphpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.hiv.primsecosyphpos <- as.numeric(sim$epi$prev.hiv.primsecosyphpos[2600, ])
round(quantile(prev.hiv.primsecosyphpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.earlysyph <- as.numeric(sim$epi$prev.earlysyph[2600, ])
round(quantile(prev.earlysyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.latesyph <- as.numeric(sim$epi$prev.latesyph[2600, ])
round(quantile(prev.latesyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.incubprim <- as.numeric(sim$epi$prev.stage.incubprim[2600, ])
round(quantile(prev.stage.incubprim, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.seco <- as.numeric(sim$epi$prev.stage.seco[2600, ])
round(quantile(prev.stage.seco, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.earlat <- as.numeric(sim$epi$prev.stage.earlat[2600, ])
round(quantile(prev.stage.earlat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.latelat <- as.numeric(sim$epi$prev.stage.latelat[2600, ])
round(quantile(prev.stage.latelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.latelatelat <- as.numeric(sim$epi$prev.stage.latelatelat[2600, ])
round(quantile(prev.stage.latelatelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.alllatelat <- as.numeric(sim$epi$prev.stage.alllatelat[2600, ])
round(quantile(prev.stage.alllatelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.tert <- as.numeric(sim$epi$prev.stage.tert[2600, ])
round(quantile(prev.stage.tert, probs = c(0.025, 0.5, 0.975)), 3)


prev.gc <- as.numeric(sim$epi$prev.gc[2600, ])
round(quantile(prev.gc, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgc <- as.numeric(sim$epi$prev.rgc[2600, ])
round(quantile(prev.rgc, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugc <- as.numeric(sim$epi$prev.ugc[2600, ])
round(quantile(prev.ugc, probs = c(0.025, 0.5, 0.975)), 3)

prev.ct <- as.numeric(sim$epi$prev.ct[2600, ])
round(quantile(prev.ct, probs = c(0.025, 0.5, 0.975)), 3)

prev.rct <- as.numeric(sim$epi$prev.rct[2600, ])
round(quantile(prev.rct, probs = c(0.025, 0.5, 0.975)), 3)

prev.uct <- as.numeric(sim$epi$prev.uct[2600, ])
round(quantile(prev.uct, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgcct <- as.numeric(sim$epi$prev.rgcct[2600, ])
round(quantile(prev.rgcct, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugcct <- as.numeric(sim$epi$prev.ugcct[2600, ])
round(quantile(prev.ugcct, probs = c(0.025, 0.5, 0.975)), 3)

# Summary of 500 sims
# Syphilis prevalence plots
par(mfrow = c(1, 1), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph.hivpos", ylab = "Prevalence", ylim = c(0, 0.15), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.primsecosyph.hivneg", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.103, col = "blue", lty = 2)
abline(h = 0.026, col = "red", lty = 2)
title("P and S Syphilis by HIV Status")
legend("topleft", c("HIV+", "HIV-"), col = c("blue", "red"), lty = c(1, 1))

plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0.00, 0.40), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "purple", qnts.col = "purple")
abline(h = 0.20, col = "blue", lty = 2)
abline(h = 0.077, col = "red", lty = 2)
abline(h = 0.2770, col = "purple", lty = 2)
title("P, S, and Early latent Syphilis Prevalence")
legend("topright", c("Primary", "Seco", "Early latent"), col = c("blue", "red", "purple"), lty = c(1, 1, 1))

plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", ylim = c(0.00, 0.60), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.44, col = "blue", lty = 2)
abline(h = 0.006, col = "red", lty = 2)
title("Late Latent and Tertiary Syphilis Prevalence")
legend("right", c("Late latent", "Tertiary"), col = c("blue", "red"), lty = c(1, 1))

# All syphilis stages
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "orange")
plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", add = TRUE, mean.col = "purple")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "black")
title("Stage-specific Syphilis Prevalence")
abline(h = 0.20, col = "blue", lty = 1)
abline(h = 0.077, col = "red", lty = 1)
abline(h = 0.2770, col = "orange", lty = 1)
abline(h = 0.44, col = "purple", lty = 1)
abline(h = 0.006, col = "black", lty = 1)
legend("topright", c("Primary/Incub", "Secondary", "Early Latent", "Late latent", "Tertiary"),
       col = c("blue", "red", "orange", "purple", "black"), lty = c(1, 1, 1, 1, 1))

# Early vs late syphilis
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.earlysyph", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.latesyph", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
title("Early and Late-stage Syphilis")
abline(h = 0.554, col = "blue", lty = 1)
abline(h = 0.446, col = "red", lty = 1)
legend("topright", c("Early", "Late"), col = c("blue", "red"), lty = c(1, 1))

# STI Prevalence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph", ylab = "Prevalence")
title("P and S Syphilis Prevalence")
abline(h = 0.046, col = "red", lty = 2)
plot(sim, y = "prev.ct", ylab = "Prevalence")
title("CT Prevalence")
plot(sim, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")

# HIV and STI Incidence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph")
abline(h = 0.9, col = "red", lty = 2)
title("Syph Incidence")



###############################
# Summary of mean simulation
# Run after burnin.process
###############################

# Syphilis prevalence plots
par(mfrow = c(1, 1), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph.hivpos", ylab = "Prevalence", ylim = c(0, 0.15), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.primsecosyph.hivneg", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.103, col = "blue", lty = 2)
abline(h = 0.026, col = "red", lty = 2)
title("P and S Syphilis by HIV Status")
legend("topleft", c("HIV+", "HIV-"), col = c("blue", "red"), lty = c(1, 1))

plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0.00, 0.40), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "purple", qnts.col = "purple")
abline(h = 0.20, col = "blue", lty = 2)
abline(h = 0.077, col = "red", lty = 2)
abline(h = 0.2770, col = "purple", lty = 2)
title("P, S, and Early latent Syphilis Prevalence")
legend("topright", c("Primary", "Seco", "Early latent"), col = c("blue", "red", "purple"), lty = c(1, 1, 1))

plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", ylim = c(0.00, 0.60), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.44, col = "blue", lty = 2)
abline(h = 0.006, col = "red", lty = 2)
title("Late Latent and Tertiary Syphilis Prevalence")
legend("right", c("Late latent", "Tertiary"), col = c("blue", "red"), lty = c(1, 1))

# All syphilis stages
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "orange")
plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", add = TRUE, mean.col = "purple")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "black")
title("Stage-specific Syphilis Prevalence")
abline(h = 0.20, col = "blue", lty = 1)
abline(h = 0.077, col = "red", lty = 1)
abline(h = 0.2770, col = "orange", lty = 1)
abline(h = 0.44, col = "purple", lty = 1)
abline(h = 0.006, col = "black", lty = 1)
legend("topright", c("Primary/Incub", "Secondary", "Early Latent", "Late latent", "Tertiary"),
       col = c("blue", "red", "orange", "purple", "black"), lty = c(1, 1, 1, 1, 1))

# Early vs late syphilis
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.earlysyph", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.latesyph", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
title("Early and Late-stage Syphilis")
abline(h = 0.554, col = "blue", lty = 1)
abline(h = 0.446, col = "red", lty = 1)
legend("topright", c("Early", "Late"), col = c("blue", "red"), lty = c(1, 1))

# STI Prevalence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph", ylab = "Prevalence")
title("P and S Syphilis Prevalence")
abline(h = 0.046, col = "red", lty = 2)
plot(sim, y = "prev.ct", ylab = "Prevalence")
title("CT Prevalence")
plot(sim, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")

# HIV and STI Incidence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph")
abline(h = 0.9, col = "red", lty = 2)
title("Syph Incidence")

df <- as.data.frame(sim)
head(df$prev.rgcct, 100)

plot(sim, y = "prev.rgcct", mean.smooth = FALSE)

df <- tail(as.data.frame(sim), 500)
sum(df$trans.main) / sum(df$incid)
sum(df$trans.casl) / sum(df$incid)
sum(df$trans.inst) / sum(df$incid)

sum(df$trans.main.gc) / sum(df$incid.gc)
sum(df$trans.casl.gc) / sum(df$incid.gc)
sum(df$trans.inst.gc) / sum(df$incid.gc)

sum(df$trans.main.ct) / sum(df$incid.ct)
sum(df$trans.casl.ct) / sum(df$incid.ct)
sum(df$trans.inst.ct) / sum(df$incid.ct)
