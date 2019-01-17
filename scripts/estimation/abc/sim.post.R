
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

prep <- readRDS(file = "data/abc.prep.rda")

out <- out_abc(wave = 5)
summary_abc(out)

coef.adj <- apply(out$param, 2, median)

coef.adj

load("fit.rda")

est[[2]]$coef.form[1:4] <- est[[2]]$coef.form[1:4] + coef.adj
dx <- netdx(est[[2]], nsims = 20, ncores = 5, nsteps = 500, verbose = FALSE)
dx

par(mfrow = c(1,1))
plot(dx)
