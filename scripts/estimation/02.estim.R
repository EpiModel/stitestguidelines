
## Paper 1 estimation file

suppressMessages(library("EpiModelHIV"))
rm(list = ls())

load("est/nwstats.20k.rda")


# 1. Main Model -----------------------------------------------------------

# Initialize network
nw.main <- base_nw_msm(st)

# Assign degree
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)

# Formulas
formation.m <- ~edges +
                nodefactor("deg.pers") +
                absdiff("sqrt.age") +
                degrange(from = 2) +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.m <- netest(nw.main,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = c(st$stats.m, 0),
                coef.diss = st$coef.diss.m,
                set.control.ergm = control.ergm(MCMLE.maxit = 250,
                                                # MPLE.max.dyad.types = 1e9,
                                                SAN.maxit = 2,
                                                SAN.nsteps.times = 2),
                verbose = TRUE)


# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw.pers <- nw.main

# Assign degree
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st)

# Formulas
formation.p <- ~edges +
                nodefactor("deg.main") +
                concurrent +
                absdiff("sqrt.age") +
                degrange(from = 3) +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.p <- netest(nw.pers,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = c(st$stats.p, 0),
                coef.diss = st$coef.diss.p,
                set.control.ergm = control.ergm(MCMLE.maxit = 250,
                                                # MPLE.max.dyad.types = 1e9,
                                                SAN.maxit = 2,
                                                SAN.nsteps.times = 2),
                verbose = TRUE)


# Fit inst model ----------------------------------------------------------

# Initialize network
nw.inst <- nw.main

# Assign degree
nw.inst <- set.vertex.attribute(nw.inst, "deg.main", nw.pers %v% "deg.main")
nw.inst <- set.vertex.attribute(nw.inst, "deg.pers", nw.main %v% "deg.pers")
table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers")

# Formulas
formation.i <- ~edges +
                nodefactor(c("deg.main", "deg.pers")) +
                nodefactor("riskg", base = 3) +
                absdiff("sqrt.age") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.i <- netest(nw.inst,
                formation = formation.i,
                target.stats = st$stats.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                set.control.ergm = control.ergm(MCMLE.maxit = 250,
                                                # MPLE.max.dyad.types = 1e9,
                                                SAN.maxit = 2,
                                                SAN.nsteps.times = 2),
                verbose = TRUE)

# Save data
est <- list(fit.m, fit.p, fit.i)
save(est, file = "est/fit.20k.rda")
