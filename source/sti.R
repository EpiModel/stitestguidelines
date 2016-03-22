
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

  rGC.sympt <- dat$attr$rGC.sympt
  uGC.sympt <- dat$attr$uGC.sympt
  rCT.sympt <- dat$attr$rCT.sympt
  uCT.sympt <- dat$attr$uCT.sympt

  ## Processes

  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive

  al <- dat$temp$al

# browser()

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
  uCT.sympt[idsInf_uct] <- rbinom(length(idsInf_uct), 1, uct.sympt.prob)

  ## Output

  # attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT

  dat$attr$rGC.sympt <- rGC.sympt
  dat$attr$uGC.sympt <- uGC.sympt
  dat$attr$rCT.sympt <- rCT.sympt
  dat$attr$uGC.sympt <- uGC.sympt

  # Summary stats
  dat$epi$incid.rgc <- length(idsInf_rgc)
  dat$epi$incid.ugc <- length(idsInf_ugc)
  dat$epi$incid.rct <- length(idsInf_rct)
  dat$epi$incid.uct <- length(idsInf_uct)


  return(dat)
}

dx.sti <- function(dat, at) {

  return(dat)
}

tx.sti <- function(dat, at) {

  return(dat)
}
