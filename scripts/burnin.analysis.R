
# stiPrEP burnin analysis

load("~/Dropbox/Projects/STIincidence/stiPrEP/scripts/burnin/data/sim.n1002.rda")

ls()

df <- as.data.frame(sim)
names(df)

plot(sim, y = c("prev.gc", "prev.rgc", "prev.ugc"), leg = TRUE)
mean(tail(df$prev.gc, 100))
mean(tail(df$prev.rgc, 100))
mean(tail(df$prev.ugc, 100))

mean(tail(df$incid.rgc))

plot(sim, y = "i.prev")
plot(sim, y = "num")
