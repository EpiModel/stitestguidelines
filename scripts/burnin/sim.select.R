
## File processing -------------------------------------------------------------

suppressMessages(library("EpiModelHIV"))
suppressMessages(library("tidyverse"))

list.files("data/")

fn <- list.files("data/", full.names = TRUE)

tdf <- data.frame(batch = NA, ir100.gc = NA, ir100.ct = NA, i.prev = NA)

for (i in 1:length(fn)) {
  load(fn[i])

  for (j in 1:sim$control$nsims) {
    df <- as.data.frame(x = sim, sim = j)
    df <- select(df, ir100.gc, ir100.ct, i.prev)
    df <- tail(df, 52)
    batch <- paste(paste(strsplit(fn[i], "[.]")[[1]][2:3], collapse = "."), j, sep = ".")
    out <- c(batch, colMeans(df))
    tdf <- rbind(tdf, out)
  }
  cat("\nFile", fn[i], "complete ...")

  if (i == length(fn)) {
    tdf <- tdf[-1, ]
    save(tdf, file = "data/tdf.rda")
  }
}


load("data/tdf.rda")
tdf_sel <- tdf[which(tdf$ir100.gc > 4.1 & tdf$ir100.gc < 4.3 &
                     tdf$ir100.ct > 6.5 & tdf$ir100.ct < 6.7 &
                     tdf$i.prev > 0.149 & tdf$i.prev < 0.151), ]
tdf_sel

load("data/sim.n3005.198.20190210.1557.rda")
ls()
s11 <- get_sims(sim, sims = 11)

df <- as.data.frame(s11)
df <- select(df, ir100.gc, ir100.ct, i.prev)
df <- tail(df, 52)


# Save as best-fitting ----------------------------------------------------
sim <- s11

save(sim, file = "est/sti.tnt.burnin.rda", compress = "xz")

