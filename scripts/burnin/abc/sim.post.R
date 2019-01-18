
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

prep <- readRDS(file = "data/abc.prep.rda")

out <- out_abc(wave = 4)
summary_abc(out)

coef.adj <- apply(out$param, 2, mean)
coef.adj

load("fit.rda")

est[[2]]$coef.form[1:4] <- est[[2]]$coef.form[1:4] + coef.adj
dx <- netdx(est[[2]], nsims = 20, ncores = 5, nsteps = 1000, verbose = FALSE)

print(dx)
plot(dx)

# EpiModel Network Diagnostics
# =======================
#   Diagnostic Method: Dynamic
# Simulations: 20
# Time Steps per Sim: 1000
#
# Formation Diagnostics
# -----------------------
#   Target Sim Mean Pct Diff Sim SD
# edges                          2022.500 2022.189    0.000 45.806
# nodefactor.deg.main.1           890.000  889.285   -0.001 35.405
# concurrent                      950.000  950.531    0.001 37.479
# absdiff.sqrt.age               1185.859 1186.838    0.001 33.736
# deg3+                                NA    0.000       NA  0.000
# offset(nodematch.role.class.I)       NA    0.000       NA  0.000
# offset(nodematch.role.class.R)       NA    0.000       NA  0.000
#
# Dissolution Diagnostics
# -----------------------
#   Target Sim Mean Pct Diff Sim SD
# Edge Duration  23.745   23.205   -0.023 22.694
# Pct Edges Diss  0.042    0.042   -0.001  0.004

setwd("../../../")
load("est/fit.rda")
est[[2]]$coef.form[1:4] <- est[[2]]$coef.form[1:4] + coef.adj
save(est, file = "est/fit.rda")
