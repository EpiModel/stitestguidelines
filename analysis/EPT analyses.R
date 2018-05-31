## EPT analyses
# Sensitivity Analyses

## Table 1 ---------------------------------------------------------

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No EPT
# Reference scenario here
load("data/sim.n8000.rda")
sim.base <- sim

incid.base <- unname(colSums(sim.base$epi$incid))
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct))

sims <- c(8000:8022)

qnt.low <- 0.25
qnt.high <- 0.75

eptcov <- rep(NA, length(sims))
eptint <- rep(NA, length(sims))
mainuptake <- rep(NA, length(sims))
persuptake <- rep(NA, length(sims))
instuptake <- rep(NA, length(sims))
mainongprov <- rep(NA, length(sims))
mainendprov <- rep(NA, length(sims))
persongprov <- rep(NA, length(sims))
persendprov <- rep(NA, length(sims))
instprov <- rep(NA, length(sims))
gctxsuccess <- rep(NA, length(sims))
cttxsuccess <- rep(NA, length(sims))

eptTx <-  rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))

gc.nnt <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))
gcct.nnt <- rep(NA, length(sims))

df <- data.frame(eptcov, eptint, mainuptake, persuptake, instuptake,
                 mainongprov, mainendprov, persongprov, persendprov, instprov,
                 gctxsuccess, cttxsuccess,

                 gc.incid, gc.pia,  gc.nnt,
                 ct.incid, ct.pia, ct.nnt,
                 gcct.incid, gcct.pia, gcct.nnt,

                 eptTx

)

for (i in seq_along(sims)) {

  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$eptcov[i] <- sim$param$ept.coverage
  df$eptint[i] <- sim$param$ept.risk.int
  df$mainuptake[i] <- sim$param$ept.uptake.partner.main
  df$persuptake[i] <- sim$param$ept.uptake.partner.pers
  df$instuptake[i] <- sim$param$ept.uptake.partner.inst
  df$mainongprov[i] <- sim$param$ept.provision.partner.main.ong
  df$mainendprov[i] <- sim$param$ept.provision.partner.main.end
  df$persongprov[i] <- sim$param$ept.provision.partner.pers.ong
  df$persendprov[i] <- sim$param$ept.provision.partner.pers.end
  df$instprov[i] <- sim$param$ept.provision.partner.inst
  df$gctxsuccess[i] <- sim$param$ept.gc.success
  df$cttxsuccess[i] <- sim$param$ept.ct.success

  # Incidence Rate over last year
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
  df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  vec.ir.gcct <- unname(colMeans(tail(sim$epi$ir100.gcct, 52)))
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")


  # PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.gc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.ct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.gcct <- unname(colSums(sim$epi$incid.gcct))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")


  # Number of partners treated in a time step
  vec.eptTx <- unname(colMeans(sim$epi$eptTx, na.rm = TRUE))
  df$eptTx[i] <- paste0(round(quantile(vec.eptTx, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(vec.eptTx, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(vec.eptTx, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")

  # Number needed to treat (need NG and CT specific provision?)
  eptdoses.gc <- unname(colSums(sim$epi$eptpartprovided_gc, na.rm = TRUE))
  eptdoses.ct <- unname(colSums(sim$epi$eptpartprovided_ct, na.rm = TRUE))
  eptdoses.gcct <- unname(colSums(sim$epi$eptpartprovided, na.rm = TRUE))

  if (is.na(mean(eptdoses.gcct))) {
    eptdoses.gc <- rep(0, 256)
    eptdoses.ct <- rep(0, 256)
    eptdoses.gcct <- rep(0, 256)
  }

  vec.gc.nnt <- (eptdoses.gc) / (incid.base.gc - incid.gc)
  vec.ct.nnt <- (eptdoses.ct) / (incid.base.ct - incid.ct)
  vec.gcct.nnt <- (eptdoses.gcct) / (incid.base.gcct - incid.gcct)

  if (is.nan(mean(vec.gcct.nnt))) {
    vec.gc.nnt <- rep(0, 256)
    vec.ct.nnt <- rep(0, 256)
    vec.gcct.nnt <- rep(0, 256)
  }

  df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                          " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                          " - ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                          ")")

  cat("*")

}

df

write.csv(df, "analysis/EPT Table 1.csv")


## EPT analyses
# Sensitivity Analyses

## Table 2 ---------------------------------------------------------

# Base - No EPT
# Reference scenario here
rm(list = ls())
load("data/sim.n8000.rda")
sim.base <- sim
sims <- c(8000:8022)

qnt.low <- 0.25
qnt.high <- 0.75

eptcov <- rep(NA, length(sims))
eptint <- rep(NA, length(sims))
mainuptake <- rep(NA, length(sims))
persuptake <- rep(NA, length(sims))
instuptake <- rep(NA, length(sims))
mainongprov <- rep(NA, length(sims))
mainendprov <- rep(NA, length(sims))
persongprov <- rep(NA, length(sims))
persendprov <- rep(NA, length(sims))
instprov <- rep(NA, length(sims))
gctxsuccess <- rep(NA, length(sims))
cttxsuccess <- rep(NA, length(sims))

eptpartelig <-  rep(NA, length(sims))
eptpartprovided <-  rep(NA, length(sims))
#eptpartuptake <-  rep(NA, length(sims))
#propindexeptElig <-  rep(NA, length(sims))
#eptprop_provided <-  rep(NA, length(sims))
eptpropprov_main <-  rep(NA, length(sims))
eptpropprov_pers <-  rep(NA, length(sims))
eptpropprov_inst <-  rep(NA, length(sims))
#eptprop_tx <-  rep(NA, length(sims))
eptuninfectedprovided <- rep(NA, length(sims))
#eptuninfecteduptake <- rep(NA, length(sims))
#eptgcinfectsti <- rep(NA, length(sims))
#eptctinfectsti <- rep(NA, length(sims))
# eptgcinfectundiaghiv <- rep(NA, length(sims))
# eptctinfectundiaghiv <- rep(NA, length(sims))
eptgcctinfectundiaghiv <- rep(NA, length(sims))
# eptgcinfecthiv <- rep(NA, length(sims))
# eptctinfecthiv <- rep(NA, length(sims))
eptgcctinfecthiv <- rep(NA, length(sims))

# eptpartprovided_gc <- rep(NA, length(sims))
# eptpartprovided_ct <- rep(NA, length(sims))
# eptpartprovided_main <- rep(NA, length(sims))
# eptpartprovided_pers <- rep(NA, length(sims))
# eptpartprovided_inst <- rep(NA, length(sims))
# eptpartuptake_main <- rep(NA, length(sims))
# eptpartuptake_pers <- rep(NA, length(sims))
# eptpartuptake_inst <- rep(NA, length(sims))
# eptpartuptake_gc <- rep(NA, length(sims))
# eptpartuptake_ct <- rep(NA, length(sims))

df <- data.frame(eptcov, eptint, mainuptake, persuptake, instuptake,
                 mainongprov, mainendprov, persongprov, persendprov, instprov,
                 gctxsuccess, cttxsuccess,

                 eptpartelig, eptpartprovided, eptpropprov_main, eptpropprov_pers, eptpropprov_inst,

                 eptuninfectedprovided, eptgcctinfecthiv, eptgcctinfectundiaghiv
                 # eptpartprovided_main, eptpartprovided_pers, eptpartprovided_inst,
                 #
                 # eptpartuptake, eptuninfecteduptake,
                 # eptpartuptake_main, eptpartuptake_pers, eptpartuptake_inst,



                 # eptpartprovided_gc, eptpartprovided_ct,
                 # eptpartuptake_gc, eptpartuptake_ct,

                 # eptgcinfectsti, eptgcinfectundiaghiv,
                 # eptctinfectsti, eptctinfectundiaghiv,

                 #propindexeptElig #eptprop_tx - this is tx success,


)

for (i in seq_along(sims)) {

  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$eptcov[i] <- sim$param$ept.coverage
  df$eptint[i] <- sim$param$ept.risk.int
  df$mainuptake[i] <- sim$param$ept.uptake.partner.main
  df$persuptake[i] <- sim$param$ept.uptake.partner.pers
  df$instuptake[i] <- sim$param$ept.uptake.partner.inst
  df$mainongprov[i] <- sim$param$ept.provision.partner.main.ong
  df$mainendprov[i] <- sim$param$ept.provision.partner.main.end
  df$persongprov[i] <- sim$param$ept.provision.partner.pers.ong
  df$persendprov[i] <- sim$param$ept.provision.partner.pers.end
  df$instprov[i] <- sim$param$ept.provision.partner.inst
  df$gctxsuccess[i] <- sim$param$ept.gc.success
  df$cttxsuccess[i] <- sim$param$ept.ct.success

  # EPT Partners Eligible
  vec.eptpartelig <- unname(colMeans(sim$epi$eptpartelig, na.rm = TRUE))
  df$eptpartelig[i] <- paste0(round(quantile(vec.eptpartelig, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.eptpartelig, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.eptpartelig, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  # EPT Partner Proportion Provided
  vec.eptpartelig <- unname(colMeans(sim$epi$eptpartelig, na.rm = TRUE))
  df$eptpartelig[i] <- paste0(round(quantile(vec.eptpartelig, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                              " (", round(quantile(vec.eptpartelig, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                              " - ", round(quantile(vec.eptpartelig, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                              ")")

  vec.mainprov <- unname(colMeans(sim$epi$eptpartprovided_main / sim$epi$eptpartelig_main, na.rm = TRUE))
  dat$epi$eptpropprov_main[i] <- paste0(round(quantile(vec.mainprov, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(vec.mainprov, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(vec.mainprov, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")
  vec.casprov <- unname(colMeans(sim$epi$eptpartprovided_pers / sim$epi$eptpartelig_main, na.rm = TRUE))
  dat$epi$eptpropprov_pers[i] <- paste0(round(quantile(vec.casprov, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(vec.casprov, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(vec.casprov, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")
  vec.instprov <- unname(colMeans(sim$epi$eptpartprovided_inst / sim$epi$eptpartelig_main, na.rm = TRUE))
  dat$epi$eptpropprov_inst[i] <- paste0(round(quantile(vec.instprov, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(vec.instprov, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(vec.instprov, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")


  # EPT Partners Provided (e.g. doses given out)
  # vec.eptpartprovided <- unname(colMeans(sim$epi$eptpartprovided, na.rm = TRUE))
  # df$eptpartprovided[i] <- paste0(round(quantile(vec.eptpartprovided, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                 " (", round(quantile(vec.eptpartprovided, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                 " - ", round(quantile(vec.eptpartprovided, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                 ")")
  # vec.eptpartprovided_gc <- unname(colMeans(sim$epi$eptpartprovided_gc, na.rm = TRUE))
  # df$eptpartprovided_gc[i] <- paste0(round(quantile(vec.eptpartprovided_gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                    " (", round(quantile(vec.eptpartprovided_gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                    " - ", round(quantile(vec.eptpartprovided_gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                    ")")
  # vec.eptpartprovided_ct <- unname(colMeans(sim$epi$eptpartprovided_ct, na.rm = TRUE))
  # df$eptpartprovided_ct[i] <- paste0(round(quantile(vec.eptpartprovided_ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                    " (", round(quantile(vec.eptpartprovided_ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                    " - ", round(quantile(vec.eptpartprovided_ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                    ")")
  # vec.eptpartprovided_main <- unname(colMeans(sim$epi$eptpartprovided_main, na.rm = TRUE))
  # df$eptpartprovided_main[i] <- paste0(round(quantile(vec.eptpartprovided_main, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                      " (", round(quantile(vec.eptpartprovided_main, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                      " - ", round(quantile(vec.eptpartprovided_main, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                      ")")
  # vec.eptpartprovided_pers <- unname(colMeans(sim$epi$eptpartprovided_pers, na.rm = TRUE))
  # df$eptpartprovided_pers[i] <- paste0(round(quantile(vec.eptpartprovided_pers, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                      " (", round(quantile(vec.eptpartprovided_pers, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                      " - ", round(quantile(vec.eptpartprovided_pers, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                      ")")
  # vec.eptpartprovided_inst <- unname(colMeans(sim$epi$eptpartprovided_inst, na.rm = TRUE))
  # df$eptpartprovided_inst[i] <- paste0(round(quantile(vec.eptpartprovided_inst, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                      " (", round(quantile(vec.eptpartprovided_inst, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                      " - ", round(quantile(vec.eptpartprovided_inst, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                      ")")

  # EPT Partners who Uptake (# of doses taken by partners)
  # vec.eptpartuptake <- unname(colMeans(sim$epi$eptpartuptake, na.rm = TRUE))
  # df$eptpartuptake[i] <- paste0(round(quantile(vec.eptpartuptake, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                               " (", round(quantile(vec.eptpartuptake, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                               " - ", round(quantile(vec.eptpartuptake, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                               ")")
  #
  # vec.eptpartuptake_gc <- unname(colMeans(sim$epi$eptpartuptake_gc, na.rm = TRUE))
  # df$eptpartuptake_gc[i] <- paste0(round(quantile(vec.eptpartuptake_gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                  " (", round(quantile(vec.eptpartuptake_gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                  " - ", round(quantile(vec.eptpartuptake_gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                  ")")
  # vec.eptpartuptake_ct <- unname(colMeans(sim$epi$eptpartuptake_ct, na.rm = TRUE))
  # df$eptpartuptake_ct[i] <- paste0(round(quantile(vec.eptpartuptake_ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                  " (", round(quantile(vec.eptpartuptake_ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                  " - ", round(quantile(vec.eptpartuptake_ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                  ")")
  # vec.eptpartuptake_main <- unname(colMeans(sim$epi$eptpartuptake_main, na.rm = TRUE))
  # df$eptpartuptake_main[i] <- paste0(round(quantile(vec.eptpartuptake_main, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                    " (", round(quantile(vec.eptpartuptake_main, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                    " - ", round(quantile(vec.eptpartuptake_main, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                    ")")
  # vec.eptpartuptake_pers <- unname(colMeans(sim$epi$eptpartuptake_pers, na.rm = TRUE))
  # df$eptpartuptake_pers[i] <- paste0(round(quantile(vec.eptpartuptake_pers, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                    " (", round(quantile(vec.eptpartuptake_pers, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                    " - ", round(quantile(vec.eptpartuptake_pers, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                    ")")
  # vec.eptpartuptake_inst <- unname(colMeans(sim$epi$eptpartuptake_inst, na.rm = TRUE))
  # df$eptpartuptake_inst[i] <- paste0(round(quantile(vec.eptpartuptake_inst, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                    " (", round(quantile(vec.eptpartuptake_inst, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                    " - ", round(quantile(vec.eptpartuptake_inst, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                    ")")

  # Proportion of Index
  # vec.propindexeptElig <- unname(colMeans(sim$epi$propindexeptElig, na.rm = TRUE))
  # df$propindexeptElig[i] <- paste0(round(quantile(vec.propindexeptElig, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                  " (", round(quantile(vec.propindexeptElig, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                  " - ", round(quantile(vec.propindexeptElig, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                  ")")

  # Proportion of Eligible Partners Provided EPT
  # vec.eptprop_provided <- unname(colMeans(sim$epi$eptpartprovided / sim$epi$eptpartelig, na.rm = TRUE))
  # df$eptprop_provided[i] <- paste0(round(quantile(vec.eptprop_provided, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                  " (", round(quantile(vec.eptprop_provided, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                  " - ", round(quantile(vec.eptprop_provided, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                  ")")

  # Proportion of eligible parners provided EPT who did not have any STI
  vec.eptuninfectedprovided <- unname(colMeans(sim$epi$eptuninfectedprovided / sim$epi$eptpartprovided, na.rm = TRUE))
  df$eptuninfectedprovided[i] <- paste0(round(quantile(vec.eptuninfectedprovided, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                        " (", round(quantile(vec.eptuninfectedprovided, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                        " - ", round(quantile(vec.eptuninfectedprovided, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                        ")")

  # Proportion of eligible parners who uptook EPT who did not have any STI
  # vec.eptuninfecteduptake <- unname(colMeans(sim$epi$eptuninfecteduptake / sim$epi$eptpartuptake, na.rm = TRUE))
  # df$eptuninfecteduptake[i] <- paste0(round(quantile(vec.eptuninfecteduptake, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                     " (", round(quantile(vec.eptuninfecteduptake, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                     " - ", round(quantile(vec.eptuninfecteduptake, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                     ")")
  # # Proportion of eligible parners provided EPT for NG who had CT, Syph, or HIV
  # vec.eptgcinfectsti <- unname(colMeans(sim$epi$eptgcinfectsti / sim$epi$eptpartuptake_gc, na.rm = TRUE))
  # df$eptgcinfectsti[i] <- paste0(round(quantile(vec.eptgcinfectsti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                " (", round(quantile(vec.eptgcinfectsti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                " - ", round(quantile(vec.eptgcinfectsti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                ")")
  # # Proportion of eligible parners provided EPT for CT who had CT, Syph, or HIV
  # vec.eptctinfectsti <- unname(colMeans(sim$epi$eptctinfectsti / sim$epi$eptpartuptake_ct, na.rm = TRUE))
  # df$eptctinfectsti[i] <- paste0(round(quantile(vec.eptctinfectsti, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                " (", round(quantile(vec.eptctinfectsti, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                " - ", round(quantile(vec.eptctinfectsti, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                ")")
  # Proportion of eligible parners provided EPT for NG who had undiagnosed HIV
  # vec.eptgcinfectundiaghiv <- unname(colMeans(sim$epi$eptgcinfectundiaghiv / sim$epi$eptpartuptake_gc, na.rm = TRUE))
  # df$eptgcinfectundiaghiv[i] <- paste0(round(quantile(vec.eptgcinfectundiaghiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                " (", round(quantile(vec.eptgcinfectundiaghiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                " - ", round(quantile(vec.eptgcinfectundiaghiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                ")")
  # # Proportion of eligible parners provided EPT for CT who had undiagnosed HIV
  # vec.eptctinfectundiaghiv <- unname(colMeans(sim$epi$eptctinfectundiaghiv / sim$epi$eptpartuptake_ct, na.rm = TRUE))
  # df$eptctinfectundiaghiv[i] <- paste0(round(quantile(vec.eptctinfectundiaghiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                " (", round(quantile(vec.eptctinfectundiaghiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                " - ", round(quantile(vec.eptctinfectundiaghiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                ")")
  # Proportion of eligible parners provided EPT for NG or CT who had undiagnosed HIV
  vec.eptgcctinfectundiaghiv <- unname(colMeans(sim$epi$eptgcctinfectundiaghiv / sim$epi$eptpartprovided, na.rm = TRUE))
  df$eptgcctinfectundiaghiv[i] <- paste0(round(quantile(vec.eptgcctinfectundiaghiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                   " (", round(quantile(vec.eptgcctinfectundiaghiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                   " - ", round(quantile(vec.eptgcctinfectundiaghiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                   ")")


  # Proportion of main parners provided EPT for NG or CT who had undiagnosed HIV
  vec.eptgcctinfectundiaghiv_main <- unname(colMeans(sim$epi$eptgcctinfectundiaghiv_main / sim$epi$eptpartprovided_main, na.rm = TRUE))
  df$eptgcctinfectundiaghiv_main[i] <- paste0(round(quantile(vec.eptgcctinfectundiaghiv_main, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                         " (", round(quantile(vec.eptgcctinfectundiaghiv_main, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                         " - ", round(quantile(vec.eptgcctinfectundiaghiv_main, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                         ")")

  # Proportion of casual parners provided EPT for NG or CT who had undiagnosed HIV
  vec.eptgcctinfectundiaghiv_pers <- unname(colMeans(sim$epi$eptgcctinfectundiaghiv_pers / sim$epi$eptpartprovided_pers, na.rm = TRUE))
  df$eptgcctinfectundiaghiv_pers[i] <- paste0(round(quantile(vec.eptgcctinfectundiaghiv_pers, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                         " (", round(quantile(vec.eptgcctinfectundiaghiv_pers, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                         " - ", round(quantile(vec.eptgcctinfectundiaghiv_pers, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                         ")")

  # Proportion of one-off parners provided EPT for NG or CT who had undiagnosed HIV
  vec.eptgcctinfectundiaghiv_inst <- unname(colMeans(sim$epi$eptgcctinfecthiv_inst / sim$epi$eptpartprovided_inst, na.rm = TRUE))
  df$eptgcctinfectundiaghiv_inst[i] <- paste0(round(quantile(vec.eptgcctinfectundiaghiv_inst, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                         " (", round(quantile(vec.eptgcctinfectundiaghiv_inst, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                         " - ", round(quantile(vec.eptgcctinfectundiaghiv_inst, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                         ")")



  # Proportion of eligible parners provided EPT for NG who had HIV
  # vec.eptgcinfecthiv <- unname(colMeans(sim$epi$eptgcinfecthiv / sim$epi$eptpartuptake_gc, na.rm = TRUE))
  # df$eptgcinfecthiv[i] <- paste0(round(quantile(vec.eptgcinfecthiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                      " (", round(quantile(vec.eptgcinfecthiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                      " - ", round(quantile(vec.eptgcinfecthiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                      ")")
  # # Proportion of eligible parners provided EPT for CT who had HIV
  # vec.eptctinfecthiv <- unname(colMeans(sim$epi$eptctinfecthiv / sim$epi$eptpartuptake_ct, na.rm = TRUE))
  # df$eptctinfecthiv[i] <- paste0(round(quantile(vec.eptctinfecthiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
  #                                      " (", round(quantile(vec.eptctinfecthiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
  #                                      " - ", round(quantile(vec.eptctinfecthiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
  #                                      ")")
  # Proportion of eligible parners provided EPT for NG or CT who had HIV
  vec.eptgcctinfecthiv <- unname(colMeans(sim$epi$eptgcctinfecthiv / sim$epi$eptpartprovided, na.rm = TRUE))
  df$eptgcctinfecthiv[i] <- paste0(round(quantile(vec.eptgcctinfecthiv, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                         " (", round(quantile(vec.eptgcctinfecthiv, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                         " - ", round(quantile(vec.eptgcctinfecthiv, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                         ")")
  # Proportion of main parners provided EPT for NG or CT who had HIV
  vec.eptgcctinfecthiv_main <- unname(colMeans(sim$epi$eptgcctinfecthiv_main / sim$epi$eptpartprovided_main, na.rm = TRUE))
  df$eptgcctinfecthiv_main[i] <- paste0(round(quantile(vec.eptgcctinfecthiv_main, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                              " (", round(quantile(vec.eptgcctinfecthiv_main, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                              " - ", round(quantile(vec.eptgcctinfecthiv_main, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                              ")")

  # Proportion of casual parners provided EPT for NG or CT who had HIV
  vec.eptgcctinfecthiv_pers <- unname(colMeans(sim$epi$eptgcctinfecthiv_pers / sim$epi$eptpartprovided_pers, na.rm = TRUE))
  df$eptgcctinfecthiv_pers[i] <- paste0(round(quantile(vec.eptgcctinfecthiv_pers, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                              " (", round(quantile(vec.eptgcctinfecthiv_pers, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                              " - ", round(quantile(vec.eptgcctinfecthiv_pers, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                              ")")

  # Proportion of one-off parners provided EPT for NG or CT who had HIV
  vec.eptgcctinfecthiv_inst <- unname(colMeans(sim$epi$eptgcctinfecthiv_inst / sim$epi$eptpartprovided_inst, na.rm = TRUE))
  df$eptgcctinfecthiv_inst[i] <- paste0(round(quantile(vec.eptgcctinfecthiv_inst, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                                              " (", round(quantile(vec.eptgcctinfecthiv_inst, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                                              " - ", round(quantile(vec.eptgcctinfecthiv_inst, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                                              ")")

  cat("*")

}

df

write.csv(df, "analysis/EPT Table 2.csv")


## Table 3 ---------------------------------------------------------

# Base - No EPT
# Reference scenario here
rm(list = ls())
load("data/sim.n8000.rda")
sim.base <- sim

incid.base <- unname(colSums(sim.base$epi$incid))
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
incid.base.gcct <- unname(colSums(sim.base$epi$incid.gcct))

# Newer way:
sims <- c(8000, 8023:8070)

qnt.low <- 0.25
qnt.high <- 0.75

eptcov <- rep(NA, length(sims))
eptint <- rep(NA, length(sims))
mainuptake <- rep(NA, length(sims))
persuptake <- rep(NA, length(sims))
instuptake <- rep(NA, length(sims))
mainongprov <- rep(NA, length(sims))
mainendprov <- rep(NA, length(sims))
persongprov <- rep(NA, length(sims))
persendprov <- rep(NA, length(sims))
instprov <- rep(NA, length(sims))
gctxsuccess <- rep(NA, length(sims))
cttxsuccess <- rep(NA, length(sims))

eptTx <-  rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.pia <- rep(NA, length(sims))
ct.incid <- rep(NA, length(sims))
ct.pia <- rep(NA, length(sims))
gcct.incid <- rep(NA, length(sims))
gcct.pia <- rep(NA, length(sims))

gc.nnt <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))
gcct.nnt <- rep(NA, length(sims))

df <- data.frame(eptcov, eptint, mainuptake, persuptake, instuptake,
                 mainongprov, mainendprov, persongprov, persendprov, instprov,
                 gctxsuccess, cttxsuccess,

                 gc.incid, gc.pia,  gc.nnt,
                 ct.incid, ct.pia, ct.nnt,
                 gcct.incid, gcct.pia, gcct.nnt,

                 eptTx

)

for (i in seq_along(sims)) {

  #fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  fn <- list.files("data/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)

  df$eptcov[i] <- sim$param$ept.coverage
  df$eptint[i] <- sim$param$ept.risk.int
  df$mainuptake[i] <- sim$param$ept.uptake.partner.main
  df$persuptake[i] <- sim$param$ept.uptake.partner.pers
  df$instuptake[i] <- sim$param$ept.uptake.partner.inst
  df$mainongprov[i] <- sim$param$ept.provision.partner.main.ong
  df$mainendprov[i] <- sim$param$ept.provision.partner.main.end
  df$persongprov[i] <- sim$param$ept.provision.partner.pers.ong
  df$persendprov[i] <- sim$param$ept.provision.partner.pers.end
  df$instprov[i] <- sim$param$ept.provision.partner.inst
  df$gctxsuccess[i] <- sim$param$ept.gc.success
  df$cttxsuccess[i] <- sim$param$ept.ct.success

  # Incidence Rate over last year
  vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
  df$gc.incid[i] <- paste0(round(quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
  df$ct.incid[i] <- paste0(round(quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")

  vec.ir.gcct <- unname(colMeans(tail(sim$epi$ir100.gcct, 52)))
  df$gcct.incid[i] <- paste0(round(quantile(vec.ir.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                             " (", round(quantile(vec.ir.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                             " - ", round(quantile(vec.ir.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                             ")")


  # PIA (Cumulative)
  incid.gc <- unname(colSums(sim$epi$incid.gc))
  vec.nia.gc <- incid.base.gc - incid.gc
  vec.pia.gc <- vec.nia.gc/incid.base.gc

  incid.ct <- unname(colSums(sim$epi$incid.ct))
  vec.nia.ct <- incid.base.ct - incid.ct
  vec.pia.ct <- vec.nia.ct/incid.base.ct

  incid.gcct <- unname(colSums(sim$epi$incid.gcct))
  vec.nia.gcct <- incid.base.gcct - incid.gcct
  vec.pia.gcct <- vec.nia.gcct/incid.base.gcct

  df$gc.pia[i] <- paste0(round(quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.pia[i] <- paste0(round(quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.pia[i] <- paste0(round(quantile(vec.pia.gcct, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.pia.gcct, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.pia.gcct, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")


  # Number of partners treated in a time step
  vec.eptTx <- unname(colMeans(sim$epi$eptTx, na.rm = TRUE))
  df$eptTx[i] <- paste0(round(quantile(vec.eptTx, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                        " (", round(quantile(vec.eptTx, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                        " - ", round(quantile(vec.eptTx, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                        ")")

  # Number needed to treat (need NG and CT specific provision?)
  eptdoses.gc <- unname(colSums(sim$epi$eptpartprovided_gc, na.rm = TRUE))
  eptdoses.ct <- unname(colSums(sim$epi$eptpartprovided_ct, na.rm = TRUE))
  eptdoses.gcct <- unname(colSums(sim$epi$eptpartprovided, na.rm = TRUE))

  if (is.na(mean(eptdoses.gcct))) {
    eptdoses.gc <- rep(0, 256)
    eptdoses.ct <- rep(0, 256)
    eptdoses.gcct <- rep(0, 256)
  }

  vec.gc.nnt <- (eptdoses.gc) / (vec.nia.gc)
  vec.ct.nnt <- (eptdoses.ct) / (vec.nia.ct)
  vec.gcct.nnt <- (eptdoses.gcct) / (vec.nia.gcct)

  if (is.nan(mean(vec.gcct.nnt))) {
    vec.gc.nnt <- rep(0, 256)
    vec.ct.nnt <- rep(0, 256)
    vec.gcct.nnt <- rep(0, 256)
  }

  df$gc.nnt[i] <- paste0(round(quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$ct.nnt[i] <- paste0(round(quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                         " (", round(quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                         " - ", round(quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                         ")")
  df$gcct.nnt[i] <- paste0(round(quantile(vec.gcct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE), 2),
                           " (", round(quantile(vec.gcct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE), 2),
                           " - ", round(quantile(vec.gcct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE), 2),
                           ")")


  cat("*")

}

df

write.csv(df, "analysis/EPT Table 3.csv")


