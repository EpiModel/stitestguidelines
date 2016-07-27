
## Process burn-in
library("EpiModelHPC")

# Examine output
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")

list.files("data/")

load("data/sim.n100.rda")
plot(sim, y = "i.prev", ylim = c(0.2, 0.3), qnts = 0.5)
abline(h = 0.26)

df <- as.data.frame(sim)
round(mean(tail(df$i.prev, 100)), 3)




# Save burn-in file for FU sims
sim <- merge_simfiles(simno = 100, indir = "data/", ftype = "max")
sim <- get_sims(sim, sims = "mean", var = "i.prev")
tail(as.data.frame(sim)$i.prev)

save(sim, file = "est/stimod.burnin.rda")
system("scp est/stimod.burnin.rda hyak:/gscratch/csde/sjenness/sti/est/")
