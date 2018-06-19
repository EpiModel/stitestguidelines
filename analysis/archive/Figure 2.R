## New version of contour plot - Figure 2
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

#source("analysis/fx.R")


# Process Data --------------------------------------------------------

#load("data/followup/sim.n3000.rda")
load("data/sim.n3000.rda")
sim.base <- sim

haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))
tests.base <- unname(colSums(sim.base$epi$hivtests.nprep))

haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))
tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests))

haz.sti <- as.numeric(colMeans(tail(sim.base$epi$ir100.sti, 52)))
ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000
incid.base.sti <- unname(colSums(sim.base$epi$incid.sti[2:521, ]))
tests.sti.base <- unname(colSums(sim.base$epi$stiasympttests))

sims <- c(3009, 3018, 3027, 3036, 3045, 3054, 3063, 3072, 3081, 3090,
          3099, 3108, 3117, 3126, 3135, 3144, 3153, 3162, 3171, 3180,
          3230:3418)

for (i in seq_along(sims)) {
  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  # PIA
  incid.gc <- unname(colSums(sim$epi$incid.gc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc
  pia.gc <- median(vec.pia.gc, na.rm = TRUE)

  incid.ct <- unname(colSums(sim$epi$incid.ct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct
  pia.ct <- median(vec.pia.ct, na.rm = TRUE)

  incid.syph <- unname(colSums(sim$epi$incid.syph))
  vec.nia.syph <- incid.base.syph - incid.syph
  vec.pia.syph <- vec.nia.syph/incid.base.syph
  pia.syph <- median(vec.pia.syph, na.rm = TRUE)

  incid.sti <- unname(colSums(sim$epi$incid.sti))
  vec.nia.sti <- incid.base.sti - incid.sti
  vec.pia.sti <- vec.nia.sti/incid.base.sti
  pia.sti <- median(vec.pia.sti, na.rm = TRUE)

  hiv.tests <- unname(colSums(sim$epi$hivtests.nprep, na.rm = TRUE))
  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))
  sti.asympt.tests <- unname(colSums(sim$epi$stiasympttests, na.rm = TRUE))

  #HIV could be HIV tests or total STI tests
  vec.nnt.hiv <- (hiv.tests - tests.base) / (incid.base - unname(colSums(sim$epi$incid)))
  vec.nnt.gc <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - unname(colSums(sim$epi$incid.gc)))
  vec.nnt.ct <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - unname(colSums(sim$epi$incid.ct)))
  vec.nnt.syph <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - unname(colSums(sim$epi$incid.syph)))
  vec.nnt.sti <- (sti.asympt.tests  - tests.sti.base) / (incid.base.sti - unname(colSums(sim$epi$incid.sti)))

  nnt.gc <- median(vec.nnt.gc, na.rm = TRUE)
  nnt.ct <- median(vec.nnt.ct, na.rm = TRUE)
  nnt.syph <- median(vec.nnt.syph, na.rm = TRUE)
  nnt.sti <- median(vec.nnt.sti, na.rm = TRUE)

  new.df <- data.frame(scenario = sims[i],
                       p1 = sim$param$stihighrisktest.ct.hivpos.coverage,
                       p2 = sim$param$partnercutoff,
                       pia.gc = pia.gc,
                       pia.ct = pia.ct,
                       pia.syph = pia.syph,
                       pia.sti = pia.sti,
                       nnt.gc = nnt.gc,
                       nnt.ct = nnt.ct,
                       nnt.syph = nnt.syph,
                       nnt.sti = nnt.sti)

  if (i == 1) {
    df <- new.df
  } else {
    df <- rbind(df, new.df)
  }

  cat("*")
}


## PIA
prev.gc.loess <- loess(pia.gc ~ p1 * p2, data = df)
prev.gc.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                              p2 = seq(1, 10, 1)))
prev.gc.fit2$PIA <- as.numeric(predict(prev.gc.loess, newdata = prev.gc.fit2))

prev.ct.loess <- loess(pia.ct ~ p1 * p2, data = df)
prev.ct.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                              p2 = seq(1, 10, 1)))
prev.ct.fit2$PIA <- as.numeric(predict(prev.ct.loess, newdata = prev.ct.fit2))

prev.syph.loess <- loess(pia.syph ~ p1 * p2, data = df)
prev.syph.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                              p2 = seq(1, 10, 1)))
prev.syph.fit2$PIA <- as.numeric(predict(prev.syph.loess, newdata = prev.syph.fit2))

prev.sti.loess <- loess(pia.sti ~ p1 * p2, data = df)
prev.sti.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                              p2 = seq(1, 10, 1)))
prev.sti.fit2$PIA <- as.numeric(predict(prev.sti.loess, newdata = prev.sti.fit2))


## NNT
prev.gc.loess.nnt <- loess(nnt.gc ~ p1 * p2, data = df)
prev.gc.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                 p2 = seq(1, 10, 1)))
prev.gc.fit2.nnt$NNT <- as.numeric(predict(prev.gc.loess.nnt, newdata = prev.gc.fit2.nnt))

prev.ct.loess.nnt <- loess(nnt.ct ~ p1 * p2, data = df)
prev.ct.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                 p2 = seq(1, 10, 1)))
prev.ct.fit2.nnt$NNT <- as.numeric(predict(prev.ct.loess.nnt, newdata = prev.ct.fit2.nnt))

prev.syph.loess.nnt <- loess(nnt.syph ~ p1 * p2, data = df)
prev.syph.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                   p2 = seq(1, 10, 1)))
prev.syph.fit2.nnt$NNT <- as.numeric(predict(prev.syph.loess.nnt, newdata = prev.syph.fit2.nnt))

prev.sti.loess.nnt <- loess(nnt.sti ~ p1 * p2, data = df)
prev.sti.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                  p2 = seq(1, 10, 1)))
prev.sti.fit2.nnt$NNT <- as.numeric(predict(prev.sti.loess.nnt, newdata = prev.sti.fit2.nnt))


# Plot --------------------------------------------------------
# tiff(filename = "analysis/Fig1b.tiff", height = 6, width = 11, units = "in", res = 250)
# plot1 <- ggplot(prev.gc.fit2, aes(p1, p2)) +
#   geom_raster(aes(fill = PIA), interpolate = TRUE) +
#   geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(title = "Percent of NG Infections Averted",
#        y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
#   # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
#   scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
#   theme(legend.position = "right")
#
# plot2 <- ggplot(prev.ct.fit2, aes(p1, p2)) +
#   geom_raster(aes(fill = PIA), interpolate = TRUE) +
#   geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(title = "Percent of CT Infections Averted",
#        y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
#   # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
#   scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
#   theme(legend.position = "right")
#
# plot3 <- ggplot(prev.syph.fit2, aes(p1, p2)) +
#   geom_raster(aes(fill = PIA), interpolate = TRUE) +
#   geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(title = "Percent of Syph Infections Averted",
#        y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
#   # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
#   scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
#   theme(legend.position = "right")
#
# plot4 <- ggplot(prev.sti.fit2, aes(p1, p2)) +
#   geom_raster(aes(fill = PIA), interpolate = TRUE) +
#   geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(title = "Percent of STI Infections Averted",
#        y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
#   # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
#   scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
#   theme(legend.position = "right")
#
# grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
#
# dev.off()


# PIA - same scale--------------------------------------------------------
tiff(filename = "analysis/Fig2PIA.tiff", height = 6, width = 11, units = "in", res = 250)

a <- rbind(prev.gc.fit2, prev.ct.fit2, prev.syph.fit2, prev.sti.fit2)
a$class[1:210] <- "Gonorrhea"
a$class[211:420] <- "Chlamydia"
a$class[421:630] <- "Syphilis"
a$class[631:840] <- "STI"

plot1 <- ggplot(a, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Infections Averted",
       y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot1)


dev.off()


# NNT - same scale--------------------------------------------------------
tiff(filename = "analysis/Fig2NNT.tiff", height = 6, width = 11, units = "in", res = 250)

b <- rbind(prev.gc.fit2.nnt, prev.ct.fit2.nnt, prev.syph.fit2.nnt, prev.sti.fit2.nnt)
b$class[1:210] <- "Gonorrhea"
b$class[211:420] <- "Chlamydia"
b$class[421:630] <- "Syphilis"
b$class[631:840] <- "STI"

plot2 <- ggplot(b, aes(p1, p2)) +
  geom_raster(aes(fill = NNT), interpolate = TRUE) +
  geom_contour(aes(z = NNT), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Number Needed to Treat",
       y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot2)


dev.off()

