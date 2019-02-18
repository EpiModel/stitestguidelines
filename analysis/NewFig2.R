## New version of contour plot - Figure 2
## SA Coverage x interval
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

# haz.gcct <- as.numeric(colMeans(tail(sim.base$epi$ir100.gcct, 52), na.rm = TRUE))
# ir.base.gcct <- unname(colMeans(sim.base$epi$ir100.gcct, na.rm = TRUE)) * 1000
# #incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct, na.rm = TRUE))
# incid.base.gcct <- unname(colSums(sim.base$epi$incid.gc, na.rm = TRUE)) + unname(colSums(sim.base$epi$incid.ct, na.rm = TRUE))
# tests.gcct.base <- unname(colSums(sim.base$epi$GCasympttests, na.rm = TRUE)) + unname(colSums(sim.base$epi$CTasympttests, na.rm = TRUE))

sims <- c(9275:9406)
for (i in seq_along(sims)) {

  fn <- list.files("data/followup/Guidelines Paper/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

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
                       p1 = sim$param$stitest.active.int,
                       p2 = sim$param$stianntest.ct.hivneg.coverage,
                       pia.gc = pia.gc,
                       pia.ct = pia.ct,
                       nnt.gc = nnt.gc,
                       nnt.ct = nnt.ct)
                        # pia.gcct = pia.gcct,
                        # nnt.gcct = nnt.gcct)
  new.df2 <- data.frame(scenario = sims[i],
                        p1 = sim$param$sti.highrisktest.int,
                        p2 = sim$param$stihighrisktest.ct.hivneg.coverage,
                        pia.gc = pia.gc,
                        pia.ct = pia.ct,
                        nnt.gc = nnt.gc,
                        nnt.ct = nnt.ct)
                        # pia.gcct = pia.gcct,
                        # nnt.gcct = nnt.gcct)


  if (i == 1) {
    df <- new.df
  } else if (i <= 45) {
    df <- rbind(df, new.df)
  } else if (i == 46) {
    df2 <- new.df2
  } else {
    df2 <- rbind(df2, new.df2)
  }

  cat("*")
}

## PIA --------------------------------------------------------
# Sexually Active Screening
prev.gc.loess <- loess(pia.gc ~ p1 * p2, data = df)
prev.gc.fit2 <- expand.grid(list(p1 = seq(26, 76, 4),
                                   p2 = seq(0.45, 0.61, 0.01)))
prev.gc.fit2$PIA <- as.numeric(predict(prev.gc.loess, newdata = prev.gc.fit2))

prev.ct.loess <- loess(pia.ct ~ p1 * p2, data = df)
prev.ct.fit2 <- expand.grid(list(p1 = seq(26, 76, 4),
                                   p2 = seq(0.45, 0.61, 0.01)))
prev.ct.fit2$PIA <- as.numeric(predict(prev.ct.loess, newdata = prev.ct.fit2))

# prev.gcct.loess <- loess(pia.gcct ~ p1 * p2, data = df)
# prev.gcct.fit2 <- expand.grid(list(p1 = seq(26, 76, 4),
#                                    p2 = seq(0.45, 0.61, 0.01)))
# prev.gcct.fit2$PIA <- as.numeric(predict(prev.gcct.loess, newdata = prev.gcct.fit2))

# HR Screening
prev.gc.loess2 <- loess(pia.gc ~ p1 * p2, data = df2)
prev.gc.fit3 <- expand.grid(list(p1 = seq(4, 52, 4),
                                   p2 = seq(0.0, 0.4, 0.05)))
prev.gc.fit3$PIA <- as.numeric(predict(prev.gc.loess2, newdata = prev.gc.fit3))

prev.ct.loess2 <- loess(pia.ct ~ p1 * p2, data = df2)
prev.ct.fit3 <- expand.grid(list(p1 = seq(4, 52, 4),
                                   p2 = seq(0.0, 0.4, 0.05)))
prev.ct.fit3$PIA <- as.numeric(predict(prev.ct.loess2, newdata = prev.ct.fit3))

# prev.gcct.loess2 <- loess(pia.gcct ~ p1 * p2, data = df2)
# prev.gcct.fit3 <- expand.grid(list(p1 = seq(4, 52, 4),
#                                    p2 = seq(0.0, 0.4, 0.05)))
# prev.gcct.fit3$PIA <- as.numeric(predict(prev.gcct.loess2, newdata = prev.gcct.fit3))

# PIA - same scale
tiff(filename = "analysis/New NGCTFig2 PIA.tiff", height = 8, width = 11, units = "in", res = 250)
par(mfrow = c(2, 2))

a <- rbind(prev.gc.fit2, prev.ct.fit2)
b <- rbind(prev.gc.fit3, prev.ct.fit3)
c <- rbind(a, b)
a$class[1:221] <- "Gonorrhea"
a$class[222:442] <- "Chlamydia"
b$class[1:117] <- "Gonorrhea"
b$class[118:234] <- "Chlamydia"
c$class <- NA
c$class[1:221] <- "Gonorrhea"
c$class[222:442] <- "Chlamydia"
c$class[443:559] <- "Gonorrhea"
c$class[560:676] <- "Chlamydia"

plot1 <- ggplot(a, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  #scale_y_discrete(labels = c("Baseline","+5%","+10%","+15%", "+20%", "+25%", "+30%", "+35%", "+40%")) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_x_discrete(labels = c("6","9","12","15", "18")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Total Infections Averted (PIA)",
       x = "Sexually Active Screening Interval (Weeks)", y = "Coverage of Sexually Active Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")


plot2 <- ggplot(b, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Total Infections Averted (PIA)",
       x = "Higher-Risk Screening Interval (Weeks)", y = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot1, plot2, nrow = 2)

dev.off()

plot3 <- ggplot(c, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Total Infections Averted (PIA)",
       x = "Screening Interval (Weeks)", y = "Coverage Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")
plot3


# NNT - same scale--------------------------------------------------------
prev.gc.loess2.nnt <- loess(nnt.gc ~ p1 * p2, data = df)
prev.gc.fit3.nnt <- expand.grid(list(p1 = seq(26, 76, 4),
                                       p2 = seq(0.45, 0.61, 0.01)))
prev.gc.fit3.nnt$NNT <- as.numeric(predict(prev.gc.loess2.nnt, newdata = prev.gc.fit3.nnt))

prev.ct.loess2.nnt <- loess(nnt.ct ~ p1 * p2, data = df)
prev.ct.fit3.nnt <- expand.grid(list(p1 = seq(26, 76, 4),
                                       p2 = seq(0.45, 0.61, 0.01)))
prev.ct.fit3.nnt$NNT <- as.numeric(predict(prev.ct.loess2.nnt, newdata = prev.ct.fit3.nnt))

# prev.gcct.loess2.nnt <- loess(nnt.gcct ~ p1 * p2, data = df)
# prev.gcct.fit3.nnt <- expand.grid(list(p1 = seq(26, 76, 4),
#                                        p2 = seq(0.45, 0.61, 0.01)))
# prev.gcct.fit3.nnt$NNT <- as.numeric(predict(prev.gcct.loess2.nnt, newdata = prev.gcct.fit3.nnt))

# HR Screening
prev.gc.loess3.nnt <- loess(nnt.gc ~ p1 * p2, data = df2)
prev.gc.fit4.nnt <- expand.grid(list(p1 = seq(4, 52, 4),
                                       p2 = seq(0.0, 0.4, 0.05)))
prev.gc.fit4.nnt$NNT <- as.numeric(predict(prev.gc.loess3.nnt, newdata = prev.gc.fit4.nnt))

prev.ct.loess3.nnt <- loess(nnt.ct ~ p1 * p2, data = df2)
prev.ct.fit4.nnt <- expand.grid(list(p1 = seq(4, 52, 4),
                                       p2 = seq(0.0, 0.4, 0.05)))
prev.ct.fit4.nnt$NNT <- as.numeric(predict(prev.ct.loess3.nnt, newdata = prev.ct.fit4.nnt))

# prev.gcct.loess3.nnt <- loess(nnt.gcct ~ p1 * p2, data = df2)
# prev.gcct.fit4.nnt <- expand.grid(list(p1 = seq(4, 52, 4),
#                                        p2 = seq(0.0, 0.4, 0.05)))
# prev.gcct.fit4.nnt$NNT <- as.numeric(predict(prev.gcct.loess3.nnt, newdata = prev.gcct.fit4.nnt))

tiff(filename = "analysis/New NGCTFig2 NNT.tiff", height = 8, width = 11, units = "in", res = 250)
par(mfrow = c(2, 2))

d <- rbind(prev.gc.fit3.nnt, prev.ct.fit3.nnt)
e <- rbind(prev.gc.fit4.nnt, prev.ct.fit4.nnt)
f <- rbind(d, e)
d$class[1:221] <- "Gonorrhea"
d$class[222:442] <- "Chlamydia"
e$class[1:117] <- "Gonorrhea"
e$class[118:234] <- "Chlamydia"
f$class <- NA
f$class[1:221] <- "Gonorrhea"
f$class[222:442] <- "Chlamydia"
f$class[443:559] <- "Gonorrhea"
f$class[560:676] <- "Chlamydia"

plot4 <- ggplot(d, aes(p1, p2)) +
  geom_raster(aes(fill = NNT), interpolate = TRUE) +
  geom_contour(aes(z = NNT), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  #scale_y_discrete(labels = c("Baseline","+5%","+10%","+15%", "+20%", "+25%", "+30%", "+35%", "+40%")) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_x_discrete(labels = c("6","9","12","15", "18")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Number Needed to Screen",
       x = "Sexually Active Screening Interval (Weeks)", y = "Coverage of Sexually Active Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")


plot5 <- ggplot(e, aes(p1, p2)) +
  geom_raster(aes(fill = NNT), interpolate = TRUE) +
  geom_contour(aes(z = NNT), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Number Needed to Screen",
       x = "Higher-Risk Screening Interval (Weeks)", y = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot4, plot5, nrow = 2)

dev.off()

plot6 <- ggplot(c, aes(p1, p2)) +
  geom_raster(aes(fill = NNT), interpolate = TRUE) +
  geom_contour(aes(z = NNT), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Number Needed to Screen",
       x = "Screening Interval (Weeks)", y = "Coverage Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")
plot6



