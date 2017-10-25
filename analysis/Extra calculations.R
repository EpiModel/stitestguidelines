sim2 <- truncate_sim(sim2, at = 5201)
sim2$epi$ir100.ct
sim2$epi$ir100.ct.tttraj1
sim2$epi$ir100.ct.tttraj2
sim2$epi$ir100.gc
sim2$epi$ir100.gc.tttraj1
sim2$epi$ir100.gc.tttraj2
sim2$epi$ir100.syph
sim2$epi$ir100.syph.tttraj1
sim2$epi$ir100.syph.tttraj2
sim2$epi$ir100.sti
sim2$epi$ir100.sti.tttraj1
sim2$epi$ir100.sti.tttraj2
sim2$epi$incid.gc
sim2$epi$incid.gc.tttraj1
sim2$epi$incid.gc.tttraj2
sim2$epi$incid.ct
sim2$epi$incid.ct.tttraj1
sim2$epi$incid.ct.tttraj2
sim2$epi$incid.syph
sim2$epi$incid.syph.tttraj1
sim2$epi$incid.syph.tttraj2
sim2$epi$incid.sti
sim2$epi$incid.sti.tttraj1
sim2$epi$incid.sti.tttraj2

# Prevalence by risk group
sim2$epi$prev.gc.tttraj1
sim2$epi$prev.gc.tttraj2
sim2$epi$prev.ct.tttraj1 #Fix
sim2$epi$prev.ct.tttraj2 #Fix
sim2$epi$prev.syph.tttraj1
sim2$epi$prev.syph.tttraj2
sim2$epi$prev.primsecosyph.tttraj1
sim2$epi$prev.primsecosyph.tttraj2
sim2$epi$prev.sti.tttraj1
sim2$epi$prev.sti.tttraj2 #Fix

# Tests by risk group
sim2$epi$rCTasympttests
sim2$epi$rCTasympttests.tttraj1
sim2$epi$rCTasympttests.tttraj2

sim2$epi$uCTasympttests
sim2$epi$uCTasympttests.tttraj1
sim2$epi$uCTasympttests.tttraj2

sim2$epi$CTasympttests
sim2$epi$CTasympttests.tttraj1
sim2$epi$CTasympttests.tttraj2

sim2$epi$rGCasympttests
sim2$epi$rGCasympttests.tttraj1
sim2$epi$rGCasympttests.tttraj2

sim2$epi$uGCasympttests
sim2$epi$uGCasympttests.tttraj1
sim2$epi$uGCasympttests.tttraj2

sim2$epi$GCasympttests
sim2$epi$GCasympttests.tttraj1
sim2$epi$GCasympttests.tttraj2

sim2$epi$syphasympttests
sim2$epi$syphasympttests.tttraj1
sim2$epi$syphasympttests.tttraj2

sim2$epi$stiasympttests
sim2$epi$stiasympttests.tttraj1
sim2$epi$stiasympttests.tttraj2

sim2$epi$rCTsympttests
sim2$epi$rCTsympttests.tttraj1
sim2$epi$rCTsympttests.tttraj2

sim2$epi$uCTsympttests
sim2$epi$uCTsympttests.tttraj1
sim2$epi$uCTsympttests.tttraj2

sim2$epi$CTsympttests
sim2$epi$CTsympttests.tttraj1
sim2$epi$CTsympttests.tttraj2

sim2$epi$rGCsympttests
sim2$epi$rGCsympttests.tttraj1
sim2$epi$rGCsympttests.tttraj2

sim2$epi$uGCsympttests
sim2$epi$uGCsympttests.tttraj1
sim2$epi$uGCsympttests.tttraj2

sim2$epi$GCsympttests
sim2$epi$GCsympttests.tttraj1
sim2$epi$GCsympttests.tttraj2

sim2$epi$syphsympttests
sim2$epi$syphsympttests.tttraj1
sim2$epi$syphsympttests.tttraj2

sim2$epi$stisympttests
sim2$epi$stisympttests.tttraj1
sim2$epi$stisympttests.tttraj2

# Treatments by risk group
sim2$epi$txGC
sim2$epi$txGC.tttraj1
sim2$epi$txGC.tttraj2

sim2$epi$txGC_asympt
sim2$epi$txGC_asympt.tttraj1
sim2$epi$txGC_asympt.tttraj2

sim2$epi$txCT
sim2$epi$txCT.tttraj1
sim2$epi$txCT.tttraj2

sim2$epi$txCT_asympt
sim2$epi$txCT_asympt.tttraj1
sim2$epi$txCT_asympt.tttraj2

sim2$epi$txsyph
sim2$epi$txsyph.tttraj1
sim2$epi$txsyph.tttraj2

sim2$epi$txsyph_asympt
sim2$epi$txsyph_asympt.tttraj1
sim2$epi$txsyph_asympt.tttraj2

sim2$epi$txearlysyph
sim2$epi$txearlysyph.tttraj1
sim2$epi$txearlysyph.tttraj2

sim2$epi$txlatesyph
sim2$epi$txlatesyph.tttraj1
sim2$epi$txlatesyph.tttraj2

sim2$epi$txSTI
sim2$epi$txSTI.tttraj1
sim2$epi$txSTI.tttraj2

sim2$epi$txSTI_asympt
sim2$epi$txSTI_asympt.tttraj1
sim2$epi$txSTI_asympt.tttraj2


# New metric
# Treatments * (persons on tx --> tx weeks)
(sum(sim2$epi$txCT) / (sim2$epi$num * sim2$epi$prev.ct)) * 52
(sum(sim2$epi$txGC) / (sim2$epi$num * sim2$epi$prev.gc)) * 52
(sum(sim2$epi$txsyph) / (sim2$epi$num * sim2$epi$prev.syph)) * 52
(sum(sim2$epi$txSTI) / (sim2$epi$num * sim2$epi$prev.sti)) * 52

(sum(sim2$epi$txCT_asympt) / (sim2$epi$num * sim2$epi$prev.ct)) * 52
(sum(sim2$epi$txGC_asympt) / (sim2$epi$num * sim2$epi$prev.gc)) * 52
(sum(sim2$epi$txsyph_asympt) / (sim2$epi$num * sim2$epi$prev.syph)) * 52
(sum(sim2$epi$txSTI_asympt) / (sim2$epi$num * sim2$epi$prev.sti)) * 52




quantile(colMeans(tail(sim$epi$prev.rgc.hivneg.only, 52)) + colMeans(tail(sim$epi$prev.ugc.hivneg.only, 52)))
quantile(colMeans(tail(sim$epi$prev.rgc.hivpos.only, 52)) + colMeans(tail(sim$epi$prev.ugc.hivpos.only, 52)))

quantile(colMeans(tail(sim$epi$prev.rct.hivneg.only, 52)) + colMeans(tail(sim$epi$prev.uct.hivneg.only, 52)))
quantile(colMeans(tail(sim$epi$prev.rct.hivpos.only, 52)) + colMeans(tail(sim$epi$prev.uct.hivpos.only, 52)))

quantile(colMeans(tail(sim$epi$prev.syph.hivneg.only, 52)))
quantile(colMeans(tail(sim$epi$prev.syph.hivpos.only, 52)))

quantile(colMeans(tail(sim$epi$prev.hivnegmultsti, 52)))
quantile(colMeans(tail(sim$epi$prev.hivposmultsti, 52)))



quantile(colMeans(tail(sim$epi$prev.rgc.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.ugc.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.rgc.hivpos.only, 52)) +
           colMeans(tail(sim$epi$prev.ugc.hivpos.only, 52)))

quantile(colMeans(tail(sim$epi$prev.rct.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.uct.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.rct.hivpos.only, 52)) +
           colMeans(tail(sim$epi$prev.uct.hivpos.only, 52)))

quantile(colMeans(tail(sim$epi$prev.syph.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.syph.hivpos.only, 52)))

quantile(colMeans(tail(sim$epi$prev.hivnegmultsti, 52)) +
           colMeans(tail(sim$epi$prev.hivposmultsti, 52)))


# HIV-
HIVneggc <- sim$epi$prev.gc.hivneg.only * (1 / (1-sim$epi$i.prev))
quantile(colMeans(tail(HIVneggc, 52)))
HIVnegct <- sim$epi$prev.ct.hivneg.only * (1 /(1-sim$epi$i.prev))
quantile(colMeans(tail(HIVnegct, 52)))
HIVnegsyph <- sim$epi$prev.syph.hivneg.only * (1 / (1-sim$epi$i.prev))
quantile(colMeans(tail(HIVnegsyph, 52)))
HIVnegmult <- sim$epi$prev.hivnegmultsti * (1 / (1-sim$epi$i.prev))
quantile(colMeans(tail(HIVnegmult, 52)))

# HIV+
HIVposgc <- sim$epi$prev.gc.hivpos.only * (1 / sim$epi$i.prev)
quantile(colMeans(tail(HIVposgc, 52)))
HIVposct <- sim$epi$prev.ct.hivpos.only * (1 / sim$epi$i.prev)
quantile(colMeans(tail(HIVposct, 52)))
HIVpossyph <- sim$epi$prev.syph.hivpos.only * (1 / sim$epi$i.prev)
quantile(colMeans(tail(HIVpossyph, 52)))
HIVposmult <- sim$epi$prev.hivposmultsti * (1 / sim$epi$i.prev)
quantile(colMeans(tail(HIVposmult, 52)))

