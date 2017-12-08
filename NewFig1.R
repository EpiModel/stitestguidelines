## New version of contour plot - Figure 1 - Syph and NG/CT
## Partner cutoff by HR coverage
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

haz.gcct <- as.numeric(colMeans(tail(sim.base$epi$ir100.gcct, 52), na.rm = TRUE))
ir.base.gcct <- unname(colMeans(sim.base$epi$ir100.gcct, na.rm = TRUE)) * 1000
#incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct, na.rm = TRUE))
incid.base.gcct <- unname(colSums(sim.base$epi$incid.gc, na.rm = TRUE)) + unname(colSums(sim.base$epi$incid.ct, na.rm = TRUE))
tests.gcct.base <- unname(colSums(sim.base$epi$GCasympttests, na.rm = TRUE)) + unname(colSums(sim.base$epi$CTasympttests, na.rm = TRUE))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52), na.rm = TRUE))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph, na.rm = TRUE)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph, na.rm = TRUE))
tests.syph.base <- unname(colSums(sim.base$epi$syphasympttests, na.rm = TRUE))


sims <- c(3009, 3018, 3027, 3036, 3045, 3054, 3063, 3072, 3081, 3090,
          3099, 3108, 3117, 3126, 3135, 3144, 3153, 3162, 3171, 3180,
          3230:3418)

for (i in seq_along(sims)) {
  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  # PIA
  #incid.gcct <- unname(colSums(sim$epi$incid.gcct, na.rm = TRUE))
  incid.gcct <- unname(colSums(sim$epi$incid.gc, na.rm = TRUE)) + unname(colSums(sim$epi$incid.ct, na.rm = TRUE))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct
  pia.gcct <- median(vec.pia.gcct, na.rm = TRUE)

  incid.syph <- unname(colSums(sim$epi$incid.syph, na.rm = TRUE))
  vec.nia.syph <- incid.base.syph - incid.syph
  vec.pia.syph <- vec.nia.syph/incid.base.syph
  pia.syph <- median(vec.pia.syph, na.rm = TRUE)

  gcct.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE)) + unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))

  vec.nnt.gcct <- (gcct.asympt.tests - tests.gcct.base) / (incid.base.gcct - incid.gcct)
  vec.nnt.syph <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - incid.syph)

  nnt.gcct <- median(vec.nnt.gcct, na.rm = TRUE)
  nnt.syph <- median(vec.nnt.syph, na.rm = TRUE)

  new.df <- data.frame(scenario = sims[i],
                       p1 = sim$param$stihighrisktest.ct.hivpos.coverage,
                       p2 = sim$param$partnercutoff,
                       pia.gcct = pia.gcct,
                       pia.syph = pia.syph,
                       nnt.gcct = nnt.gcct,
                       nnt.syph = nnt.syph)

  if (i == 1) {
    df <- new.df
  } else {
    df <- rbind(df, new.df)
  }

  cat("*")
}


## PIA
prev.gcct.loess <- loess(pia.gcct ~ p1 * p2, data = df)
prev.gcct.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                 p2 = seq(1, 10, 1)))
prev.gcct.fit2$PIA <- as.numeric(predict(prev.gcct.loess, newdata = prev.gcct.fit2))

prev.syph.loess <- loess(pia.syph ~ p1 * p2, data = df)
prev.syph.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                   p2 = seq(1, 10, 1)))
prev.syph.fit2$PIA <- as.numeric(predict(prev.syph.loess, newdata = prev.syph.fit2))


## NNT
prev.gcct.loess.nnt <- loess(nnt.gcct ~ p1 * p2, data = df)
prev.gcct.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                     p2 = seq(1, 10, 1)))
prev.gcct.fit2.nnt$NNT <- as.numeric(predict(prev.gcct.loess.nnt, newdata = prev.gcct.fit2.nnt))

prev.syph.loess.nnt <- loess(nnt.syph ~ p1 * p2, data = df)
prev.syph.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                       p2 = seq(1, 10, 1)))
prev.syph.fit2.nnt$NNT <- as.numeric(predict(prev.syph.loess.nnt, newdata = prev.syph.fit2.nnt))

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
tiff(filename = "analysis/New Fig1 PIA.tiff", height = 8, width = 11, units = "in", res = 250)

a <- rbind(prev.gcct.fit2, prev.syph.fit2)
a$class[1:210] <- "Gonorrhea/Chlamydia"
a$class[211:420] <- "Syphilis"

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
tiff(filename = "analysis/New Fig1 NNT.tiff", height = 8, width = 11, units = "in", res = 250)

b <- rbind(prev.gcct.fit2.nnt, prev.syph.fit2.nnt)
b$class[1:210] <- "Gonorrhea/Chlamydia"
b$class[211:420] <- "Syphilis"

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

