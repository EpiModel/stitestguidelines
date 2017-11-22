## New version of contour plot - Figure 3
## Coverage x interval
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

sims <- c(3419:3423, 3189:3193, 3424:3458, 3459:3463, 3194:3198, 3464:3513)


for (i in seq_along(sims)) {
  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  incid.gcct <- unname(colSums(sim$epi$incid.gcct))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct
  pia.gcct <- median(vec.pia.gcct, na.rm = TRUE)

  incid.syph <- unname(colSums(sim$epi$incid.syph))
  vec.nia.syph <- incid.base.syph - incid.syph
  vec.pia.syph <- vec.nia.syph/incid.base.syph
  pia.syph <- median(vec.pia.syph, na.rm = TRUE)

  gcct.asympt.tests <- unname(colSums(sim$epi$GCasympttests, na.rm = TRUE)) + unname(colSums(sim$epi$CTasympttests, na.rm = TRUE))
  syph.asympt.tests <- unname(colSums(sim$epi$syphasympttests, na.rm = TRUE))

  vec.nnt.gcct <- (gcct.asympt.tests - tests.gcct.base) / (incid.base.gcct - unname(colSums(sim$epi$incid.gcct)))
  vec.nnt.syph <- (syph.asympt.tests  - tests.syph.base) / (incid.base.syph - unname(colSums(sim$epi$incid.syph)))

  nnt.gcct <- median(vec.nnt.gcct, na.rm = TRUE)
  nnt.syph <- median(vec.nnt.syph, na.rm = TRUE)

  new.df <- data.frame(scenario = sims[i],
                       p1 = sim$param$stitest.active.int,
                       p2 = sim$param$stianntest.ct.hivneg.coverage,
                       pia.gcct = pia.gcct,
                       pia.syph = pia.syph,
                       nnt.gcct = nnt.gcct,
                       nnt.syph = nnt.syph)
  new.df2 <- data.frame(scenario = sims[i],
                        p1 = sim$param$sti.highrisktest.int,
                        p2 = sim$param$stihighrisktest.ct.hivneg.coverage,
                        pia.gcct = pia.gcct,
                        pia.syph = pia.syph,
                        nnt.gcct = nnt.gcct,
                        nnt.syph = nnt.syph)


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

#rownames(df) <-
#rownames(df2) <-


## PIA
# Sexually Active Screening
prev.gcct.loess <- loess(pia.gcct ~ p1 * p2, data = df)
prev.gcct.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                   p2 = seq(1, 10, 1)))
prev.gcct.fit2$PIA <- as.numeric(predict(prev.gcct.loess, newdata = prev.gcct.fit2))

prev.syph.loess <- loess(pia.syph ~ p1 * p2, data = df)
prev.syph.fit2 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                   p2 = seq(1, 10, 1)))
prev.syph.fit2$PIA <- as.numeric(predict(prev.syph.loess, newdata = prev.syph.fit2))

# prev.sti.fit2 <- expand.grid(list(p1 = c(26, 39, 52, 64, 77),
#                                   p2 = seq(0.44, 0.616, 0.022)))

# HR Screening
prev.gcct.loess2 <- loess(pia.gcct ~ p1 * p2, data = df2)
prev.gcct.fit3 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                   p2 = seq(1, 10, 1)))
prev.gcct.fit3$PIA <- as.numeric(predict(prev.gcct.loess2, newdata = prev.gcct.fit2))

prev.syph.loess2 <- loess(pia.syph ~ p1 * p2, data = df2)
prev.syph.fit3 <- expand.grid(list(p1 = seq(0.0, 1, 0.05),
                                   p2 = seq(1, 10, 1)))
prev.syph.fit3$PIA <- as.numeric(predict(prev.syph.loess2, newdata = prev.syph.fit2))

# prev.sti.fit3 <- expand.grid(list(p1 = c(4, 13, 26, 39, 52),
#                                   p2 = seq(0.0, 1.0, 0.05)))



# PIA - same scale--------------------------------------------------------
tiff(filename = "analysis/New Fig2 PIA.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(2, 1))

a <- rbind(prev.gcct.fit2, prev.syph.fit2)
b <- rbind(prev.gcct.fit3, prev.syph.fit3)
c <- rbind(a, b)
a$class[1:45] <- "STI"
b$class[1:105] <- "STI"
c$class <- NA
c$class[1:45] <- "SA"
c$class[46:150] <- "HR"

plot1 <- ggplot(a, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  #scale_y_discrete(labels = c("Baseline","+5%","+10%","+15%", "+20%", "+25%", "+30%", "+35%", "+40%")) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_x_discrete(labels = c("6","9","12","15", "18")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Infections Averted",
       x = "Screening Interval (Weeks)", y = "Coverage of Sexually Active Screening") +
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
  labs(title = "Percent of Infections Averted",
       x = "Screening Interval (Weeks)", y = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot1, plot2, ncol = 2)

dev.off()

plot3 <- ggplot(c, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  facet_wrap(~class, scales = 'fixed', ncol = 2, nrow = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Infections Averted",
       x = "Screening Interval (Weeks)", y = "Coverage Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "Spectral", direction = -1) +
  theme(legend.position = "right")
plot3




