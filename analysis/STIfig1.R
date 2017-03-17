# Figure 1 ----------------------------------------------------------------

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise", n = 9, type = "continuous")

# Baseline 
load("data/sim.n3000.rda")
sim.base <- truncate_sim(sim, at = 2600)
mn.base <- as.data.frame(sim.base)
ir.gc.base <- (sum(mn.base$incid.gc)/sum((1 - mn.base$prev.gc) * mn.base$num)) * 52 * 1e5
ir.ct.base <- (sum(mn.base$incid.ct)/sum((1 - mn.base$prev.ct) * mn.base$num)) * 52 * 1e5
ir.syph.base <- (sum(mn.base$incid.syph)/sum((1 - mn.base$prev.syph) * mn.base$num)) * 52 * 1e5
incid.gc.base <- sum(mn.base$incid.gc)
incid.ct.base <- sum(mn.base$incid.ct)
incid.syph.base <- sum(mn.base$incid.syph)

# Sims 3142-3152 for comparisons to all - 6 month interval
# 3142 - STI
# 3143 - recent partners
# 3144 - new partners
# 3145 - partner who has multiple partners
# 3146 - partner with a STI
# 3147 - any CAI in a non-main
# 3148 - any CAI
# 3149 - recent or new partners
# 3150 - sti, recent, or new partners
# 3151 - CAI in non-main or any CAI
# 3152 - partner with multiple partners or with a STI

sims <- c(3142:3152)
df.gc.pia <- data.frame(rep(NA, 256))
df.gc.nnt <- data.frame(rep(NA, 256))
df.ct.pia <- data.frame(rep(NA, 256))
df.ct.nnt <- data.frame(rep(NA, 256))
df.syph.pia <- data.frame(rep(NA, 256))
df.syph.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {
    
    load(list.files("data/", pattern = as.character(sims[i]), full.names = TRUE))
    sim <- truncate_sim(sim, at = 2600)
    mn <- as.data.frame(sim)
    ir.gc <- (colSums(sim$epi$incid.gc, na.rm = TRUE)) /
        sum((1 - mn$prev.gc)  * mn$num) * 52 * 1e5
    ir.ct <- (colSums(sim$epi$incid.ct, na.rm = TRUE)) /
        sum((1 - mn$prev.ct)  * mn$num) * 52 * 1e5
    ir.syph <- (colSums(sim$epi$incid.syph, na.rm = TRUE)) /
        sum((1 - mn$prev.syph)  * mn$num) * 52 * 1e5
    
    vec.gc.nia <- round(ir.gc.base - unname(ir.gc), 1)
    df.gc.pia[, i] <- vec.gc.nia / ir.gc.base
    
    vec.ct.nia <- round(ir.ct.base - unname(ir.ct), 1)
    df.ct.pia[, i] <- vec.ct.nia / ir.ct.base
    
    vec.syph.nia <- round(ir.syph.base - unname(ir.syph), 1)
    df.syph.pia[, i] <- vec.syph.nia / ir.syph.base
    

    gc.tests <- unname(tail(sim$epi$totalGCasympttests, 1))
    ct.tests <- unname(tail(sim$epi$totalCTasympttests, 1))
    syph.tests <- unname(tail(sim$epi$totalsyphasympttests, 1))

    df.gc.nnt[, i] <- gc.tests/(incid.gc.base - unname(colSums(sim$epi$incid.gc)))
    df.ct.nnt[, i] <- ct.tests/(incid.ct.base - unname(colSums(sim$epi$incid.ct)))
    df.syph.nnt[, i] <- syph.tests/(incid.syph.base - unname(colSums(sim$epi$incid.syph)))
    
}
names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "J1", "J2", "J3", "J4")

head(df.gc.pia)
head(df.ct.pia)
head(df.syph.pia)
head(df.gc.nnt)
head(df.ct.nnt)
head(df.syph.nnt)

pal <- wes_palette("Zissou")[c(1, 5)]
# pdf(file = "analysis/P1Fig2.pdf", height = 6, width = 12, pointsize = 16)
tiff(filename = "analysis/Fig1.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))

# Left Panel: PIA
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "PIA by Behavioral Indication",
        xlab = "Behavioral Indication", ylab = "Percent NG Infections Averted")
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "PIA by Behavioral Indication",
        xlab = "Behavioral Indication", ylab = "Percent CT Infections Averted")
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "PIA by Behavioral Indication",
        xlab = "Behavioral Indication", ylab = "Percent Syph Infections Averted")

# Right Panel: NNT
boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "NNT by Behavioral Indication",
        xlab = "Behavioral Indication", ylab = "NG Number Needed to Treat")
boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "NNT by Behavioral Indication",
        xlab = "Behavioral Indication", ylab = "CT Number Needed to Treat")
boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 6), rep(pal[2], 3)),
        main = "NNT by Behavioral Indication",
        xlab = "Behavioral Indication", ylab = "Syph Number Needed to Treat")

dev.off()
