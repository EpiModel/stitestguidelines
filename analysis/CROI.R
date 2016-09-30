
# CROI Table ----------------------------------------------------------

sim.base <- merge_simfiles(1000)
sim.comp <- merge_simfiles(1001)

epi_stats(sim.base, at = 520, qnt.low = 0.025, qnt.high = 0.975)
epi_stats(sim.base, sim.comp, at = 520, qnt.low = 0.025, qnt.high = 0.975)


# Loading baseline simulations (No PrEP)
load("data/sim.n1300.rda")
sim.base <- truncate_sim(sim, at = 2600)

# 20% coverage
load("data/sim.n1303.rda")
sim.comp <- truncate_sim(sim, at = 2600)

sim.comp$param$prep.coverage
sim.comp$param$prep.sti.screen.int
sim.comp$param$rcomp.prob

df.base <- as.data.frame(sim.base)
df.base$txCT

df.comp <- as.data.frame(sim.comp)
df.comp$txCT

load("data/sim.n1404.rda")
sim.base <- truncate_sim(sim, at = 2600)

load("data/sim.n1405.rda")
sim.comp <- truncate_sim(sim, at = 2600)

df.base <- as.data.frame(sim.base)
df.comp <- as.data.frame(sim.comp)
tail(df.base$prop.CT.asympt.tx)
tail(sim.base$epi$prop.CT.asympt.tx)

mean(tail(df.base$prop.rGC.tx, 104))
mean(tail(df.base$prop.rCT.tx, 104))

mean(tail(df.comp$prop.rGC.tx, 104))
mean(tail(df.comp$prop.rCT.tx, 104))

mean(tail(df.comp$prop.rGC.tx, 104))/mean(tail(df.base$prop.rGC.tx, 104))
mean(tail(df.comp$prop.rCT.tx, 104))/mean(tail(df.base$prop.rCT.tx, 104))

mean(tail(df.comp$prop.CT.asympt.tx, 104))
mean(tail(df.comp$prop.GC.asympt.tx, 104))
