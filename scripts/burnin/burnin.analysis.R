
# stiPrEP burnin analysis

rm(list=ls())

list.files("data/")
unlink("data/sim.*")
system("scp hyak:/gscratch/csde/sjenness/sti/data/sim.n300.rda scripts/burnin/data/")
load("scripts/burnin/data/sim.n300.rda")

ls()

df <- tail(as.data.frame(sim), 250)
mean(df$prev.rgcct)
mean(df$prev.ugcct)
mean(df$ir100.gc)
mean(df$ir100.ct)
mean(df$i.prev)

par(mfrow = c(1,1))
plot(sim, y = "i.prev", ylim = c(0.20, 0.30))
plot(sim, y = "prev.rgcct", sim.lines = TRUE)
plot(sim, y = "prev.ugcct")
plot(sim, y = "ir100.gc")
plot(sim, y = "ir100.ct")

df <- as.data.frame(sim)
head(df$prev.rgcct, 100)

plot(sim, y = "prev.rgcct", mean.smooth = FALSE)

df <- tail(as.data.frame(sim), 500)
sum(df$trans.main) / sum(df$incid)
sum(df$trans.casl) / sum(df$incid)
sum(df$trans.inst) / sum(df$incid)
