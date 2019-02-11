
library("methods")
suppressMessages(library("EpiABC"))
suppressMessages(library("EpiModel"))

out <- get_posterior(wave = 15, input = "scripts/estimation/abc/data/")
summary(out)

plot(out, type = "stats")
plot(out, type = "param")

coef.adj <- colMeans(out$param)
coef.adj

load("est/fit.20k.rda")

est[[2]]$coef.form[1:4] <- est[[2]]$coef.form[1:4] + coef.adj
dx <- netdx(est[[2]], nsims = 25, ncores = 25, nsteps = 1000, verbose = TRUE)

print(dx)
plot(dx)

# EpiModel Network Diagnostics
# =======================
#   Diagnostic Method: Dynamic
# Simulations: 25
# Time Steps per Sim: 1000
#
# Formation Diagnostics
# -----------------------
#   Target Sim Mean Pct Diff Sim SD
# edges                          4045.000 4050.702    0.001 68.757
# nodefactor.deg.main.1          1780.000 1778.343   -0.001 47.536
# concurrent                     1900.000 1901.081    0.001 55.357
# absdiff.sqrt.age               2371.718 2374.343    0.001 49.575
# deg3+                                NA    0.000       NA  0.000
# offset(nodematch.role.class.I)       NA    0.000       NA  0.000
# offset(nodematch.role.class.R)       NA    0.000       NA  0.000
#
# Dissolution Diagnostics
# -----------------------
#   Target Sim Mean Pct Diff Sim SD
# Edge Duration  23.745   23.204   -0.023 22.698
# Pct Edges Diss  0.042    0.042    0.000  0.003

save(est, file = "est/fit.20k.rda")
