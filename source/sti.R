
## STI Modules for stiPrEP Model

## TODO

sti_trans <- function(dat, at) {

  ## Parameters
  rgc.tprob <- dat$param$rgc.tprob
  ugc.tprob <- dat$param$ugc.tprob
  rct.tprob <- dat$param$rct.tprob
  uct.tprob <- dat$param$uct.tprob

  rgc.sympt.prob <- dat$param$rgc.sympt.prob
  ugc.sympt.prob <- dat$param$ugc.sympt.prob
  rct.sympt.prob <- dat$param$rct.sympt.prob
  uct.sympt.prob <- dat$param$uct.sympt.prob

  sti.cond.rr <- dat$param$sti.cond.rr

  ## Attributes
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  rGC.infTime <- dat$attr$rGC.infTime
  uGC.infTime <- dat$attr$uGC.infTime
  rCT.infTime <- dat$attr$rCT.infTime
  uCT.infTime <- dat$attr$uCT.infTime

  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt


  # set disease status to 0 for new births
  newBirths <- which(dat$attr$arrival.time == at)
  rGC[newBirths] <- 0
  uGC[newBirths] <- 0
  rCT[newBirths] <- 0
  uCT[newBirths] <- 0

  ## Processes

  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive

  al <- dat$temp$al

  # rectal GC infection
  # requires urethral GC in infected partner, infected insertive, and no rGC in sus partner
  p1Inf_rgc <- which(uGC[al[, "p1"]] == 1 & uGC[al[, "p2"]] == 0 &
                     rGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rgc <- which(uGC[al[, "p1"]] == 0 & uGC[al[, "p2"]] == 1 &
                     rGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rgc <- c(p1Inf_rgc, p2Inf_rgc)

  uai_rgc <- al[, "uai"][allActs_rgc]
  tprob_rgc <- rep(rgc.tprob, length(allActs_rgc))
  tprob_rgc[uai_rgc == 0] <- tprob_rgc[uai_rgc == 0] * sti.cond.rr

  trans_rgc <- rbinom(length(allActs_rgc), 1, tprob_rgc)

  transAL_rgc <- al[allActs_rgc[trans_rgc == 1], ]
  idsInf_rgc <- unique(ifelse(uGC[transAL_rgc[, "p1"]] == 1,
                              transAL_rgc[, "p2"], transAL_rgc[, "p1"]))

  rGC[idsInf_rgc] <- 1
  rGC.infTime[idsInf_rgc] <- at
  rGC.sympt[idsInf_rgc] <- rbinom(length(idsInf_rgc), 1, rgc.sympt.prob)


  # urethral GC infection
  # requires rectal GC in infected partner, infected receptive, and no urthethralGC in sus partner
  p1Inf_ugc <- which(rGC[al[, "p1"]] == 1 & rGC[al[, "p2"]] == 0 &
                     uGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_ugc <- which(rGC[al[, "p1"]] == 0 & rGC[al[, "p2"]] == 1 &
                     uGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_ugc <- c(p1Inf_ugc, p2Inf_ugc)

  uai_ugc <- al[, "uai"][allActs_ugc]
  tprob_ugc <- rep(ugc.tprob, length(allActs_ugc))
  tprob_ugc[uai_ugc == 0] <- tprob_ugc[uai_ugc == 0] * sti.cond.rr

  trans_ugc <- rbinom(length(allActs_ugc), 1, tprob_ugc)

  transAL_ugc <- al[allActs_ugc[trans_ugc == 1], ]
  idsInf_ugc <- unique(ifelse(uGC[transAL_ugc[, "p1"]] == 1,
                              transAL_ugc[, "p2"], transAL_ugc[, "p1"]))

  uGC[idsInf_ugc] <- 1
  uGC.infTime[idsInf_ugc] <- at
  uGC.sympt[idsInf_ugc] <- rbinom(length(idsInf_ugc), 1, ugc.sympt.prob)


  # rectal CT infection
  # requires urethral CT in infected partner, infected insertive, and no rCT in sus partner
  p1Inf_rct <- which(uCT[al[, "p1"]] == 1 & uCT[al[, "p2"]] == 0 &
                     rCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf_rct <- which(uCT[al[, "p1"]] == 0 & uCT[al[, "p2"]] == 1 &
                     rCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs_rct <- c(p1Inf_rct, p2Inf_rct)

  uai_rct <- al[, "uai"][allActs_rct]
  tprob_rct <- rep(rct.tprob, length(allActs_rct))
  tprob_rct[uai_rct == 0] <- tprob_rct[uai_rct == 0] * sti.cond.rr

  trans_rct <- rbinom(length(allActs_rct), 1, tprob_rct)

  transAL_rct <- al[allActs_rct[trans_rct == 1], ]
  idsInf_rct <- unique(ifelse(uCT[transAL_rct[, "p1"]] == 1,
                              transAL_rct[, "p2"], transAL_rct[, "p1"]))

  rCT[idsInf_rct] <- 1
  rCT.infTime[idsInf_rct] <- at
  rCT.sympt[idsInf_rct] <- rbinom(length(idsInf_rct), 1, rct.sympt.prob)

  # urethral CT infection
  # requires rectal CT in infected partner, infected receptive, and no urthethralCT in sus partner
  p1Inf_uct <- which(rCT[al[, "p1"]] == 1 & rCT[al[, "p2"]] == 0 &
                     uCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf_uct <- which(rCT[al[, "p1"]] == 0 & rCT[al[, "p2"]] == 1 &
                     uCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs_uct <- c(p1Inf_uct, p2Inf_uct)

  uai_uct <- al[, "uai"][allActs_uct]
  tprob_uct <- rep(uct.tprob, length(allActs_uct))
  tprob_uct[uai_uct == 0] <- tprob_uct[uai_uct == 0] * sti.cond.rr

  trans_uct <- rbinom(length(allActs_uct), 1, tprob_uct)

  transAL_uct <- al[allActs_uct[trans_uct == 1], ]
  idsInf_uct <- unique(ifelse(uCT[transAL_uct[, "p1"]] == 1,
                              transAL_uct[, "p2"], transAL_uct[, "p1"]))

  uCT[idsInf_uct] <- 1
  uCT.infTime[idsInf_uct] <- at
  uCT.sympt[idsInf_uct] <- rbinom(length(idsInf_uct), 1, uct.sympt.prob)

  ## Output

  # attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT

  dat$attr$rGC.infTime <- rGC.infTime
  dat$attr$uGC.infTime <- uGC.infTime
  dat$attr$rCT.infTime <- rCT.infTime
  dat$attr$uCT.infTime <- uCT.infTime

  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$rCT.sympt <- rCT.sympt
  dat$attr$uCT.sympt <- uCT.sympt

  # Summary stats
  dat$epi$incid.rgc <- length(idsInf_rgc)
  dat$epi$incid.ugc <- length(idsInf_ugc)
  dat$epi$incid.rct <- length(idsInf_rct)
  dat$epi$incid.uct <- length(idsInf_uct)

  stopifnot(all(!is.na(rGC.infTime[rGC == 1])),
            all(!is.na(rGC.sympt[rGC == 1])),
            all(!is.na(uGC.infTime[uGC == 1])),
            all(!is.na(uGC.sympt[uGC == 1])),
            all(!is.na(rCT.infTime[rCT == 1])),
            all(!is.na(rCT.sympt[rCT == 1])),
            all(!is.na(uCT.infTime[uCT == 1])),
            all(!is.na(uCT.sympt[uCT == 1])))

  return(dat)
}


sti_recov <- function(dat, at) {

  # parameters
  rgc.dur.asympt <- dat$param$rgc.dur.asympt
  ugc.dur.asympt <- dat$param$ugc.dur.asympt
  gc.dur.tx <- dat$param$gc.dur.tx
  gc.dur.ntx <- dat$param$gc.dur.ntx

  rct.dur.asympt <- dat$param$rct.dur.asympt
  uct.dur.asympt <- dat$param$uct.dur.asympt
  ct.dur.tx <- dat$param$ct.dur.tx
  ct.dur.ntx <- dat$param$ct.dur.ntx

  # attributes


  # GC recovery
  idsRGC_asympt <- which(dat$attr$rGC == 1 & dat$attr$rGC.infTime < at &
                         dat$attr$rGC.sympt == 0)
  idsUGC_asympt <- which(dat$attr$uGC == 1 & dat$attr$uGC.infTime < at &
                         dat$attr$uGC.sympt == 0)
  idsGC_tx <- which(dat$attr$rGC == 1 & dat$attr$rGC.infTime < at &
                    dat$attr$rGC.sympt == 1 & dat$attr$rGC.tx == 1)
  idsGC_ntx <- which(dat$attr$rGC == 1 & dat$attr$rGC.infTime < at &
                     dat$attr$rGC.sympt == 1 & dat$attr$rGC.tx == 0)

  recovRGC_asympt <- idsRGC_asympt[which(rbinom(length(idsRGC_asympt), 1, 1/rgc.dur.asympt) == 1)]
  recovUGC_asympt <- idsUGC_asympt[which(rbinom(length(idsUGC_asympt), 1, 1/ugc.dur.asympt) == 1)]
  recovGC_tx <- idsGC_tx[which(rbinom(length(idsGC_tx), 1, 1/gc.dur.tx) == 1)]
  recovGC_ntx <- idsGC_ntx[which(rbinom(length(idsGC_ntx), 1, 1/gc.dur.ntx) == 1)]

  recovRGC_sympt <- intersect(which(dat$attr$rGC == 1), c(recovGC_tx, recovGC_ntx))
  recovUGC_sympt <- intersect(which(dat$attr$uGC == 1), c(recovGC_tx, recovGC_ntx))

  recovRGC <- c(recovRGC_asympt, recovRGC_sympt)
  recovUGC <- c(recovUGC_asympt, recovUGC_sympt)

  dat$attr$rGC[recovRGC] <- 0
  dat$attr$rGC.sympt[recovRGC] <- NA
  dat$attr$rGC.infTime[recovRGC] <- NA
  dat$attr$rGC.tx[recovRGC] <- NA

  dat$attr$uGC[recovUGC] <- 0
  dat$attr$uGC.sympt[recovUGC] <- NA
  dat$attr$uGC.infTime[recovUGC] <- NA
  dat$attr$uGC.tx[recovUGC] <- NA


  # CT recovery
  idsRCT_asympt <- which(dat$attr$rCT == 1 & dat$attr$rCT.infTime < at &
                         dat$attr$rCT.sympt == 0)
  idsUCT_asympt <- which(dat$attr$uCT == 1 & dat$attr$uCT.infTime < at &
                         dat$attr$uCT.sympt == 0)
  idsCT_tx <- which(dat$attr$rCT == 1 & dat$attr$rCT.infTime < at &
                    dat$attr$rCT.sympt == 1 & dat$attr$rCT.tx == 1)
  idsCT_ntx <- which(dat$attr$rCT == 1 & dat$attr$rCT.infTime < at &
                     dat$attr$rCT.sympt == 1 & dat$attr$rCT.tx == 0)

  recovRCT_asympt <- idsRCT_asympt[which(rbinom(length(idsRCT_asympt), 1, 1/rct.dur.asympt) == 1)]
  recovUCT_asympt <- idsUCT_asympt[which(rbinom(length(idsUCT_asympt), 1, 1/uct.dur.asympt) == 1)]
  recovCT_tx <- idsCT_tx[which(rbinom(length(idsCT_tx), 1, 1/ct.dur.tx) == 1)]
  recovCT_ntx <- idsCT_ntx[which(rbinom(length(idsCT_ntx), 1, 1/ct.dur.ntx) == 1)]

  recovRCT_sympt <- intersect(which(dat$attr$rCT == 1), c(recovCT_tx, recovCT_ntx))
  recovUCT_sympt <- intersect(which(dat$attr$uCT == 1), c(recovCT_tx, recovCT_ntx))

  recovRCT <- c(recovRCT_asympt, recovRCT_sympt)
  recovUCT <- c(recovUCT_asympt, recovUCT_sympt)

  dat$attr$rCT[recovRCT] <- 0
  dat$attr$rCT.sympt[recovRCT] <- NA
  dat$attr$rCT.infTime[recovRCT] <- NA
  dat$attr$rCT.tx[recovRCT] <- NA

  dat$attr$uCT[recovUCT] <- 0
  dat$attr$uCT.sympt[recovUCT] <- NA
  dat$attr$uCT.infTime[recovUCT] <- NA
  dat$attr$uCT.tx[recovUCT] <- NA

  return(dat)
}


sti_tx <- function(dat, at) {

  # params
  gc.prob.tx <- dat$param$gc.prob.tx
  ct.prob.tx <- dat$param$ct.prob.tx

  # gc treatment
  idsRGC_tx <- which(dat$attr$rGC == 1 & dat$attr$rGC.infTime < at & dat$attr$rGC.sympt == 1 & is.na(dat$attr$rGC.tx))
  idsUGC_tx <- which(dat$attr$uGC == 1 & dat$attr$uGC.infTime < at & dat$attr$uGC.sympt == 1 & is.na(dat$attr$uGC.tx))
  idsGC_tx <- c(idsRGC_tx, idsUGC_tx)

  txGC <- idsGC_tx[which(rbinom(length(idsGC_tx), 1, gc.prob.tx) == 1)]
  txRGC <- intersect(idsRGC_tx, txGC)
  txUGC <- intersect(idsUGC_tx, txGC)

  # ct treatment
  idsRCT_tx <- which(dat$attr$rCT == 1 & dat$attr$rCT.infTime < at & dat$attr$rCT.sympt == 1 & is.na(dat$attr$rCT.tx))
  idsUCT_tx <- which(dat$attr$uCT == 1 & dat$attr$uCT.infTime < at & dat$attr$uCT.sympt == 1 & is.na(dat$attr$uCT.tx))
  idsCT_tx <- c(idsRCT_tx, idsUCT_tx)

  txCT <- idsCT_tx[which(rbinom(length(idsCT_tx), 1, ct.prob.tx) == 1)]
  txRCT <- intersect(idsRCT_tx, txCT)
  txUCT <- intersect(idsUCT_tx, txCT)

  # update attr
  dat$attr$rGC.tx[idsRGC_tx] <- 0
  dat$attr$rGC.tx[txRGC] <- 1

  dat$attr$uGC.tx[idsUGC_tx] <- 0
  dat$attr$uGC.tx[txUGC] <- 1

  dat$attr$rCT.tx[idsRCT_tx] <- 0
  dat$attr$rCT.tx[txRCT] <- 1

  dat$attr$uCT.tx[idsUCT_tx] <- 0
  dat$attr$uCT.tx[txUCT] <- 1

  return(dat)
}
