
# stiPrEP burnin analysis

# rm(list=ls())
# 
# list.files("data/")
# unlink("data/sim.*")
# system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n300.rda scripts/burnin/data/")
# load("scripts/burnin/data/sim.n300.rda")
# 
# ls()

## Prevalence of HIV, GC, CT, Syphilis
i.prev <- as.numeric(mean(sim$epi$i.prev[2500:2600, ]))
round(quantile(i.prev, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph <- as.numeric(mean(sim$epi$prev.syph[2500:2600, ]))
round(quantile(prev.syph, probs = c(0.025, 0.5, 0.975)), 3)

prev.gc <- as.numeric(mean(sim$epi$prev.gc[2500:2600, ]))
round(quantile(prev.gc, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgc <- as.numeric(mean(sim$epi$prev.rgc[2500:2600, ]))
round(quantile(prev.rgc, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugc <- as.numeric(mean(sim$epi$prev.ugc[2500:2600, ]))
round(quantile(prev.ugc, probs = c(0.025, 0.5, 0.975)), 3)

prev.rct <- as.numeric(mean(sim$epi$prev.rct[2500:2600, ]))
round(quantile(prev.rct, probs = c(0.025, 0.5, 0.975)), 3)

prev.uct <- as.numeric(mean(sim$epi$prev.uct[2500:2600, ]))
round(quantile(prev.uct, probs = c(0.025, 0.5, 0.975)), 3)

prev.ct <- as.numeric(mean(sim$epi$prev.ct[2500:2600, ]))
round(quantile(prev.ct, probs = c(0.025, 0.5, 0.975)), 3)


## Incidence

ir100.gc <- as.numeric(mean(sim$epi$ir100.gc[2500:2600, ]))
round(quantile(ir100.gc, probs = c(0.025, 0.5, 0.975)), 3)

ir100.ct <- as.numeric(mean(sim$epi$ir100.ct[2500:2600, ]))
round(quantile(ir100.ct, probs = c(0.025, 0.5, 0.975)), 3)

ir100 <- as.numeric(mean(sim$epi$ir100[2500:2600, ]))
round(quantile(ir100, probs = c(0.025, 0.5, 0.975)), 3)

ir100.syph <- as.numeric(mean(sim$epi$ir100.syph[2500:2600, ]))
round(quantile(ir100.syph, probs = c(0.025, 0.5, 0.975)), 3)


# Syphilis stage prevalence and HIV-syphilis overlap

prev.hiv.syphpos <- as.numeric(mean(sim$epi$prev.hiv.syphpos[2500:2600, ]))
round(quantile(prev.hiv.syphpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph.hivpos <- as.numeric(mean(sim$epi$prev.syph.hivpos[2500:2600, ]))
round(quantile(prev.syph.hivpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph.hivneg <- as.numeric(mean(sim$epi$prev.syph.hivneg[2500:2600, ]))
round(quantile(prev.syph.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.incubprim <- as.numeric(mean(sim$epi$prev.stage.incubprim[2500:2600, ]))
round(quantile(prev.stage.incubprim, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.seco <- as.numeric(mean(sim$epi$prev.stage.seco[2500:2600, ]))
round(quantile(prev.stage.seco, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.earlat <- as.numeric(mean(sim$epi$prev.stage.earlat[2500:2600, ]))
round(quantile(prev.stage.earlat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.latelat <- as.numeric(mean(sim$epi$prev.stage.latelat[2500:2600, ]))
round(quantile(prev.stage.latelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.latelatelat <- as.numeric(mean(sim$epi$prev.stage.latelatelat[2500:2600, ]))
round(quantile(prev.stage.latelatelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.alllatelat <- as.numeric(mean(sim$epi$prev.stage.alllatelat[2500:2600, ]))
round(quantile(prev.stage.alllatelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.tert <- as.numeric(mean(sim$epi$prev.stage.tert[2500:2600, ]))
round(quantile(prev.stage.tert, probs = c(0.025, 0.5, 0.975)), 3)

prev.earlysyph <- as.numeric(mean(sim$epi$prev.earlysyph[2500:2600, ]))
round(quantile(prev.earlysyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.latesyph <- as.numeric(mean(sim$epi$prev.latesyph[2500:2600, ]))
round(quantile(prev.latesyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgcct <- as.numeric(sim$epi$prev.rgcct[2600, ])
round(quantile(prev.rgcct, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugcct <- as.numeric(sim$epi$prev.ugcct[2600, ])
round(quantile(prev.ugcct, probs = c(0.025, 0.5, 0.975)), 3)

# Summary of 500 sims
par(mfrow = c(3,2), oma=c(0,0,2,0))
plot(sim, y = "i.prev", ylim = c(0.23, 0.29))
title("HIV Prevalence")
plot(sim, y = "prev.rgcct", ylim = c(0.12, 0.17))
title("Rectal STI Prevalence")
plot(sim, y = "prev.ugcct", ylim = c(0.06, 0.11))
title("Urethral STI Prevalence")
plot(sim, y = "ir100.gc", ylim = c(18, 28))
title("Gonococcal Incidence Rate")
plot(sim, y = "ir100.ct", ylim = c(20, 30))
title("Chlamydial Incidence Rate")
plot(sim, y = "ir100", ylim = c(3, 5))
title("HIV Incidence Rate")
title(main="Statistics for 500 Simulations", outer=TRUE)


# Summary of mean simulation
# Run after burnin.process
par(mfrow = c(3,2), oma=c(0,0,2,0))
plot(sim, y = "i.prev", ylim = c(0.23, 0.29))
title("HIV Prevalence")
plot(sim, y = "prev.rgcct", ylim = c(0.12, 0.17))
title("Rectal STI Prevalence")
plot(sim, y = "prev.ugcct", ylim = c(0.06, 0.11))
title("Urethral STI Prevalence")
plot(sim, y = "ir100.gc", ylim = c(18, 28))
title("Gonococcal Incidence Rate")
plot(sim, y = "ir100.ct", ylim = c(20, 30))
title("Chlamydial Incidence Rate")
plot(sim, y = "ir100", ylim = c(3, 5))
title("HIV Incidence Rate")
title(main="Statistics for Best-fitting Simulation", outer=TRUE)

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
