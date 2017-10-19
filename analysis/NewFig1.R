## New version of contour plot
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")

#source("analysis/fx.R")


# Process Data --------------------------------------------------------

load("data/followup/sim.n3000.rda")
sim.base <- sim
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
ir.base.sti <- unname(colMeans(sim.base$epi$ir100.sti)) * 1000

sims <- c(3009, 3018, 3027, 3036, 3045, 3054, 3063, 3072, 3081, 3090,
          3099, 3108, 3117, 3126, 3135, 3144, 3153, 3162, 3171, 3180,
          3230:3418)

for (i in seq_along(sims)) {
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  # PIA
  ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
  vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
  vec.pia.gc <- ifelse(vec.nia.gc > 0, vec.nia.gc/ir.base.gc, 0)
  vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
  pia.gc <- median(vec.pia.gc, na.rm = TRUE)

  ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
  vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
  vec.pia.ct <- ifelse(vec.nia.ct > 0, vec.nia.ct/ir.base.ct, 0)
  vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
  pia.ct <- median(vec.pia.ct, na.rm = TRUE)

  ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
  vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
  vec.pia.syph <- ifelse(vec.nia.syph > 0, vec.nia.syph/ir.base.syph, 0)
  vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
  pia.syph <- median(vec.pia.syph, na.rm = TRUE)

  ir.comp.sti <- unname(colMeans(sim$epi$ir100.sti)) * 1000
  vec.nia.sti <- round(ir.base.sti - ir.comp.sti, 1)
  vec.pia.sti <- ifelse(vec.nia.sti > 0, vec.nia.sti/ir.base.sti, 0)
  vec.pia.sti <- vec.pia.sti[vec.pia.sti > -Inf]
  pia.sti <- median(vec.pia.sti, na.rm = TRUE)

  new.df <- data.frame(scenario = sims[i],
                       p1 = sim$param$stihighrisktest.ct.hivpos.coverage,
                       p2 = sim$param$partnercut,
                       pia.gc = pia.gc,
                       pia.ct = pia.ct,
                       pia.syph = pia.syph,
                       pia.sti = pia.sti)

  if (i == 1) {
    df <- new.df
  } else {
    df <- rbind(df, new.df)
  }

  cat("*")
}

nrow(df)

table(df$p1)
table(df$p2)
table(df$p1, df$p2)

#df$p2 <- rep(seq(0.0, 1, 0.05), times = 16)

table(df$p1, df$p2)

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

tiff(filename = "analysis/Fig1b.tiff", height = 6, width = 11, units = "in", res = 250)

plot1 <- ggplot(prev.gc.fit2, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of NG Infections Averted",
       y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
  theme(legend.position = "right")

plot2 <- ggplot(prev.ct.fit2, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of CT Infections Averted",
       y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
  theme(legend.position = "right")

plot3 <- ggplot(prev.syph.fit2, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of Syph Infections Averted",
       y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
  theme(legend.position = "right")

plot4 <- ggplot(prev.sti.fit2, aes(p1, p2)) +
  geom_raster(aes(fill = PIA), interpolate = TRUE) +
  geom_contour(aes(z = PIA), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Percent of STI Infections Averted",
       y = "Partner Number Cutoff", x = "Coverage of Higher-Risk Screening") +
  # scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1) +
  scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
  theme(legend.position = "right")

grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)

dev.off()
