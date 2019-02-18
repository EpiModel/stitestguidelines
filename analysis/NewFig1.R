## New version of contour plot - Figure 1
## Partner cutoff by HR coverage
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

# Process Data --------------------------------------------------------

load("data/followup/Guidelines Paper/sim.n9000.rda")
sim.base <- sim

haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52), na.rm = TRUE))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc, na.rm = TRUE)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.rgc, na.rm = TRUE)) + unname(colSums(sim.base$epi$incid.ugc, na.rm = TRUE))
tests.gc.base <- unname(colSums(sim.base$epi$GCasympttests, na.rm = TRUE))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52), na.rm = TRUE))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct, na.rm = TRUE)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.rct, na.rm = TRUE)) + unname(colSums(sim.base$epi$incid.uct, na.rm = TRUE))
tests.ct.base <- unname(colSums(sim.base$epi$CTasympttests, na.rm = TRUE))

sims <- c(9000, 9009:9028, 9046:9234)

for (i in seq_along(sims)) {

  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  # PIA
  incid.gc <- unname(colSums(sim$epi$incid.rgc, na.rm = TRUE)) + unname(colSums(sim$epi$incid.ugc, na.rm = TRUE))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc
  pia.gc <- median(vec.pia.gc, na.rm = TRUE)

  incid.ct <- unname(colSums(sim$epi$incid.rct, na.rm = TRUE)) + unname(colSums(sim$epi$incid.uct, na.rm = TRUE))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct
  pia.ct <- median(vec.pia.ct, na.rm = TRUE)

  # incid.gcct <- unname(colSums(sim$epi$incid.gc, na.rm = TRUE)) + unname(colSums(sim$epi$incid.ct, na.rm = TRUE))
  # vec.nia.gcct <- incid.base.gcct - incid.gcct
  # vec.pia.gcct <- vec.nia.gcct/incid.base.gcct
  # pia.gcct <- median(vec.pia.gcct, na.rm = TRUE)

  gc.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE))
  ct.asympt.tests <- unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  # gcct.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE)) + unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))

  vec.nnt.gc <- (gc.asympt.tests - tests.gc.base) / (incid.base.gc - incid.gc)
  vec.nnt.ct <- (ct.asympt.tests - tests.ct.base) / (incid.base.ct - incid.ct)
  # vec.nnt.gcct <- (gcct.asympt.tests - tests.gcct.base) / (incid.base.gcct - incid.gcct)

  nnt.gc <- median(vec.nnt.gc, na.rm = TRUE)
  nnt.ct <- median(vec.nnt.ct, na.rm = TRUE)
  # nnt.gcct <- median(vec.nnt.gcct, na.rm = TRUE)

  new.df <- data.frame(scenario = sims[i],
                       p1 = sim$param$stihighrisktest.ct.hivpos.coverage,
                       p2 = sim$param$partnercutoff,
                       pia.gc = pia.gc,
                       pia.ct = pia.ct,
                       nnt.gc = nnt.gc,
                       nnt.ct = nnt.ct)
                       # pia.gcct = pia.gcct,
                       # nnt.gcct = nnt.gcct)

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

## NNT
prev.gc.loess.nnt <- loess(nnt.gc ~ p1 * p2, data = df)
prev.gc.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                       p2 = seq(1, 10, 1)))
prev.gc.fit2.nnt$NNT <- as.numeric(predict(prev.gc.loess.nnt, newdata = prev.gc.fit2.nnt))

prev.ct.loess.nnt <- loess(nnt.ct ~ p1 * p2, data = df)
prev.ct.fit2.nnt <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                       p2 = seq(1, 10, 1)))
prev.ct.fit2.nnt$NNT <- as.numeric(predict(prev.ct.loess.nnt, newdata = prev.ct.fit2.nnt))


# PIA - same scale--------------------------------------------------------
tiff(filename = "analysis/New NGCTFig1 PIA.tiff", height = 8, width = 11, units = "in", res = 250)

a <- rbind(prev.gc.fit2, prev.ct.fit2)
a$class[1:210] <- "Gonorrhea"
a$class[211:420] <- "Chlamydia"

plot1 <- ggplot(a, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Total Infections Averted (PIA)",
       y = "Partner Number Threshold for Higher-Risk Screening",
       x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot1)


dev.off()


# NNT - same scale--------------------------------------------------------
tiff(filename = "analysis/New NGCTFig1 NNT.tiff", height = 8, width = 11, units = "in", res = 250)

b <- rbind(prev.gc.fit2.nnt, prev.ct.fit2.nnt)
b$class[1:210] <- "Gonorrhea"
b$class[211:420] <- "Chlamydia"

plot2 <- ggplot(b, aes(p1, p2)) +
  geom_raster(aes(fill = NNT), interpolate = TRUE) +
  geom_contour(aes(z = NNT), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Number Needed to Treat",
       y = "Partner Number Threshold for Higher-Risk Screening",
       x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot2)


dev.off()

