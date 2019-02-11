
library("EpiModelHIV")

list.files("data/")
# unlink(list.files("data/", full.names = TRUE))

# sim <- merge_simfiles(simno = 1000)
# sim
load("data/sim.n2000.rda")

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(1,3))
plot(sim, y = "ir100.ct", ylim = c(0, 12), mean.smooth = FALSE)
abline(h = 6.6)
mean(colMeans(tail(sim$epi$ir100.ct, 52)))
plot(sim, y = "ir100.gc", ylim = c(0, 12), mean.smooth = FALSE)
abline(h = 4.2)
mean(colMeans(tail(sim$epi$ir100.gc, 52)))
plot(sim, y = "i.prev", mean.smooth = FALSE)
abline(h = 0.15)
mean(colMeans(tail(sim$epi$i.prev, 52)))


years <- c(40, 45, 50)
gc <- sapply(years, function(x) mean(rowMeans(sim$epi$ir100.gc[(x*52):(x*52+51), ])))
ct <- sapply(years, function(x) mean(rowMeans(sim$epi$ir100.ct[(x*52):(x*52+51), ])))
hiv <- sapply(years, function(x) mean(rowMeans(sim$epi$i.prev[(x*52):(x*52+51), ])))

out <- c(gc, ct, hiv)



setdiff(vars, names(sim$epi))

sim$epi <- sim$epi[names(sim$epi) %in% vars]
names(sim$epi)

vars <- c(
  # HIV
  "incid", "hivtests.nprep", "i.prev",

  # GC
  "incid.gc", "incid.gc.tttraj1", "incid.gc.tttraj2",
  "ir100.gc", "ir100.gc.tttraj1", "ir100.gc.tttraj2",
  "prev.gc", "prev.gc.tttraj1", "prev.gc.tttraj2",
  "GCasympttests", "GCsympttests",
  "GCasympttests.tttraj1", "GCasympttests.tttraj2",
  "GCsympttests.tttraj1", "GCsympttests.tttraj2",
  "txGC", "txGC.tttraj1", "txGC.tttraj2",
  "txGC_asympt",
  "tt.traj.gc1", "tt.traj.gc2",

  # CT
  "incid.ct", "incid.ct.tttraj1", "incid.ct.tttraj2",
  "ir100.ct", "ir100.ct.tttraj1", "ir100.ct.tttraj2",
  "prev.ct", "prev.ct.tttraj1", "prev.ct.tttraj2",
  "CTasympttests", "CTsympttests",
  "CTasympttests.tttraj1", "CTasympttests.tttraj2",
  "CTsympttests.tttraj1", "CTsympttests.tttraj2",
  "txCT", "txCT.tttraj1", "txCT.tttraj2",
  "txCT_asympt",
  "tt.traj.ct1", "tt.traj.ct2",

  # Combined
  "incid.gcct", "incid.gcct.tttraj1", "incid.gcct.tttraj2",
  "ir100.gcct", "ir100.gcct.tttraj1", "ir100.gcct.tttraj2",
  "prev.sti", "prev.sti.tttraj1", "prev.sti.tttraj2",
  "stiasympttests", "stisympttests",
  "stiasympttests.tttraj1", "stiasympttests.tttraj2",
  "stisympttests.tttraj1", "stisympttests.tttraj2",
  "txSTI", "txSTI.tttraj1", "txSTI.tttraj2",
  "txSTI_asympt",
  "tt.traj.sti1", "tt.traj.sti2",

  # Other
  "num")
