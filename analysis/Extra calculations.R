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

# Old PAF-related calcs ------------------------------------
# NG only
quantile(colMeans(tail(sim2$epi$prev.rgc.hivneg.only, 52)) + colMeans(tail(sim2$epi$prev.ugc.hivneg.only, 52)))
quantile(colMeans(tail(sim2$epi$prev.rgc.hivpos.only, 52)) + colMeans(tail(sim2$epi$prev.ugc.hivpos.only, 52)))

# CT only
quantile(colMeans(tail(sim2$epi$prev.rct.hivneg.only, 52)) + colMeans(tail(sim2$epi$prev.uct.hivneg.only, 52)))
quantile(colMeans(tail(sim2$epi$prev.rct.hivpos.only, 52)) + colMeans(tail(sim2$epi$prev.uct.hivpos.only, 52)))

# Syph only
quantile(colMeans(tail(sim2$epi$prev.syph.hivneg.only, 52)))
quantile(colMeans(tail(sim2$epi$prev.syph.hivpos.only, 52)))

# Any syphilis
quantile(colMeans(tail(sim2$epi$prev.primsecosyph.hivneg, 52)))
quantile(colMeans(tail(sim2$epi$prev.primsecosyph.hivpos, 52)))

# Multiple STI
quantile(colMeans(tail(sim2$epi$prev.hivnegmultsti, 52)))
quantile(colMeans(tail(sim2$epi$prev.hivposmultsti, 52)))

# Discordant Edges
quantile(colMeans(tail(sim2$epi$prop.edges.negneg, 52)))
quantile(colMeans(tail(sim2$epi$prop.edges.negpos, 52)))
quantile(colMeans(tail(sim2$epi$prop.edges.pospos, 52)))

# Syph IR by HIV status
quantile(colMeans(tail(sim2$epi$ir100.syph.hivneg, 52)))
quantile(colMeans(tail(sim2$epi$ir100.syph.hivpos, 52)))

# Prevalence of HIV among IPS syphilis positive
quantile(colMeans(tail(sim2$epi$prev.hiv.primsecosyphpos, 52)))

# Rectal and urethral STI
quantile(colMeans(tail(sim2$epi$prev.rgcct, 52)))
quantile(colMeans(tail(sim2$epi$prev.ugcct, 52)))


quantile(colMeans(tail(sim$epi$prev.rgc.hivneg.only, 52)) + colMeans(tail(sim$epi$prev.ugc.hivneg.only, 52)))



# PAF-related calcs ------------------------------------
# NG only
quantile(colMeans(tail(sim$epi$prev.rgc.hivneg, 52)) + colMeans(tail(sim$epi$prev.ugc.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.rgc.hivpos, 52)) + colMeans(tail(sim$epi$prev.ugc.hivpos, 52)))

# CT only
quantile(colMeans(tail(sim$epi$prev.rct.hivneg, 52)) + colMeans(tail(sim$epi$prev.uct.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.rct.hivpos, 52)) + colMeans(tail(sim$epi$prev.uct.hivpos, 52)))

# Any syphilis
quantile(colMeans(tail(sim$epi$prev.primsecosyph.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.primsecosyph.hivpos, 52)))

# Multiple STI
quantile(colMeans(tail(sim$epi$prev.hivnegmultsti, 52)))
quantile(colMeans(tail(sim$epi$prev.hivposmultsti, 52)))

# Number of edges
quantile(colMeans(tail(sim$epi$prop.edges.negneg, 52)))
quantile(colMeans(tail(sim$epi$prop.edges.negpos, 52)))
quantile(colMeans(tail(sim$epi$prop.edges.pospos, 52)))

# Number of acts by partnership serostatuses
quantile(colMeans(tail(sim$epi$num.acts.negneg, 52)))
quantile(colMeans(tail(sim$epi$num.acts.negpos, 52)))
quantile(colMeans(tail(sim$epi$num.acts.pospos, 52)))

# Proportion of UAI acts by partnership serostatuses
quantile(colMeans(tail(sim$epi$prop.uai.negneg, 52)))
quantile(colMeans(tail(sim$epi$prop.uai.negpos, 52)))
quantile(colMeans(tail(sim$epi$prop.uai.pospos, 52)))

# Proportion of acts by partnership serostatuses
quantile(colMeans(tail(sim$epi$prop.acts.negneg, 52)))
quantile(colMeans(tail(sim$epi$prop.acts.negpos, 52)))
quantile(colMeans(tail(sim$epi$prop.acts.pospos, 52)))

# Syph IR by HIV status
quantile(colMeans(tail(sim$epi$ir100.syph.hivneg, 52)))
quantile(colMeans(tail(sim$epi$ir100.syph.hivpos, 52)))

# Prevalence of HIV among IPS syphilis positive
quantile(colMeans(tail(sim$epi$prev.hiv.primsecosyphpos, 52)))
quantile(colMeans(tail(sim$epi$prev.hiv.primsecosyphneg, 52)))

# Prevalence of diagnosed HIV among IPS syphilis positive
quantile(colMeans(tail(sim$epi$prev.dxhiv.dxipssyph, 52)))

# Rectal and urethral STI
quantile(colMeans(tail(sim$epi$prev.rgc.hivneg, 52)) + colMeans(tail(sim$epi$prev.rct.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.rgc.hivpos, 52)) + colMeans(tail(sim$epi$prev.rct.hivpos, 52)))
quantile(colMeans(tail(sim$epi$prev.ugc.hivneg, 52)) + colMeans(tail(sim$epi$prev.uct.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.ugc.hivpos, 52)) + colMeans(tail(sim$epi$prev.uct.hivpos, 52)))

# Total GC prevalence
quantile(colMeans(tail(sim$epi$prev.rgc.hivneg, 52)) +
           colMeans(tail(sim$epi$prev.ugc.hivneg, 52)) +
           colMeans(tail(sim$epi$prev.rgc.hivpos, 52)) +
           colMeans(tail(sim$epi$prev.ugc.hivpos, 52)))

# Total CT prevalence
quantile(colMeans(tail(sim$epi$prev.rct.hivneg, 52)) +
           colMeans(tail(sim$epi$prev.uct.hivneg, 52)) +
           colMeans(tail(sim$epi$prev.rct.hivpos, 52)) +
           colMeans(tail(sim$epi$prev.uct.hivpos, 52)))

# Total syph prevalence
quantile(colMeans(tail(sim$epi$prev.syph.hivneg, 52)) +
           colMeans(tail(sim$epi$prev.syph.hivpos, 52)))

# Total multiple STI prevalence
quantile(colMeans(tail(sim$epi$prev.hivnegmultsti, 52)) +
           colMeans(tail(sim$epi$prev.hivposmultsti, 52)))



#NG
median(colSums(sim$epi$GCasympttests)) + median(colSums(sim$epi$GCsympttests)) #290,518.5
median(colSums(sim$epi$txGC)) #867.5

j = 1
disc.rate = 0.03
(sum(sim$epi$GCasympttests[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$GCasympttests[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$GCasympttests[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$GCasympttests[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$GCasympttests[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$GCasympttests[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$GCasympttests[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$GCasympttests[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$GCasympttests[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$GCasympttests[470:521, j]) * ((1 - disc.rate)^9)) # 244759.2
(sum(sim$epi$GCsympttests[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$GCsympttests[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$GCsympttests[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$GCsympttests[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$GCsympttests[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$GCsympttests[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$GCsympttests[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$GCsympttests[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$GCsympttests[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$GCsympttests[470:521, j]) * ((1 - disc.rate)^9)) # 916.47
# Total discounted tests = 245,675.7

(sum(sim$epi$txGC[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$txGC[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$txGC[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$txGC[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$txGC[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$txGC[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$txGC[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$txGC[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$txGC[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$txGC[470:521, j]) * ((1 - disc.rate)^9))
# Total discounted treatments = 1333.371


#CT
median(colSums(sim$epi$CTasympttests)) + median(colSums(sim$epi$CTsympttests)) #271,745.5
median(colSums(sim$epi$txCT)) #985

j = 1
disc.rate = 0.03
(sum(sim$epi$CTasympttests[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$CTasympttests[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$CTasympttests[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$CTasympttests[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$CTasympttests[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$CTasympttests[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$CTasympttests[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$CTasympttests[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$CTasympttests[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$CTasympttests[470:521, j]) * ((1 - disc.rate)^9)) # 226764.7
(sum(sim$epi$CTsympttests[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$CTsympttests[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$CTsympttests[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$CTsympttests[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$CTsympttests[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$CTsympttests[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$CTsympttests[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$CTsympttests[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$CTsympttests[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$CTsympttests[470:521, j]) * ((1 - disc.rate)^9)) # 358.1929
# Total discounted tests = 227122.9

(sum(sim$epi$txCT[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$txCT[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$txCT[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$txCT[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$txCT[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$txCT[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$txCT[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$txCT[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$txCT[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$txCT[470:521, j]) * ((1 - disc.rate)^9))
# Total discounted treatments = 1039.149


#Syph
median(colSums(sim$epi$syphasympttests)) + median(colSums(sim$epi$syphsympttests)) # 196,192
median(colSums(sim$epi$txearlysyph)) # 1564
median(colSums(sim$epi$txlatesyph)) #109.5

j = 1
disc.rate = 0.03
(sum(sim$epi$syphasympttests[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$syphasympttests[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$syphasympttests[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$syphasympttests[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$syphasympttests[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$syphasympttests[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$syphasympttests[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$syphasympttests[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$syphasympttests[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$syphasympttests[470:521, j]) * ((1 - disc.rate)^9)) # 165538.6
(sum(sim$epi$syphsympttests[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$syphsympttests[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$syphsympttests[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$syphsympttests[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$syphsympttests[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$syphsympttests[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$syphsympttests[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$syphsympttests[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$syphsympttests[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$syphsympttests[470:521, j]) * ((1 - disc.rate)^9)) # 2407.39
# Total discounted tests = 167,946

(sum(sim$epi$txearlysyph[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$txearlysyph[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$txearlysyph[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$txearlysyph[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$txearlysyph[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$txearlysyph[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$txearlysyph[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$txearlysyph[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$txearlysyph[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$txearlysyph[470:521, j]) * ((1 - disc.rate)^9))
# Total discounted treatments = 2957.019

(sum(sim$epi$txlatesyph[2:53, j]) * ((1 - disc.rate)^0) +
    sum(sim$epi$txlatesyph[54:105, j]) * ((1 - disc.rate)^1) +
    sum(sim$epi$txlatesyph[106:157, j]) * ((1 - disc.rate)^2) +
    sum(sim$epi$txlatesyph[158:209, j]) * ((1 - disc.rate)^3) +
    sum(sim$epi$txlatesyph[210:261, j]) * ((1 - disc.rate)^4) +
    sum(sim$epi$txlatesyph[262:313, j]) * ((1 - disc.rate)^5) +
    sum(sim$epi$txlatesyph[314:365, j]) * ((1 - disc.rate)^6) +
    sum(sim$epi$txlatesyph[366:417, j]) * ((1 - disc.rate)^7) +
    sum(sim$epi$txlatesyph[418:469, j]) * ((1 - disc.rate)^8) +
    sum(sim$epi$txlatesyph[470:521, j]) * ((1 - disc.rate)^9))
# Total discounted treatments = 123.808
