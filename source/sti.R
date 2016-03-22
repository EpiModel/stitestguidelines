
## STI Modules for stiPrEP Model

## TODO

sti_trans <- function(dat, at) {

  ## Parameters
  gc.rt.patp <- dat$param$gc.rt.patp
  gc.ur.patp <- dat$param$gc.ur.patp
  ct.rt.patp <- dat$param$ct.rt.patp
  ct.ur.patp <- dat$param$ct.ur.patp

  sti.cond.rr <- dat$param$sti.cond.rr

  ## Attributes
  rectalGC <- dat$attr$rectalGC
  urethralGC <- dat$attr$urethralGC
  rectalCT <- dat$attr$rectalCT
  urethralCT <- dat$attr$urethralCT

  ## Processes

  # ins = 0 : p2 is insertive
  # ins = 1 : p1 is insertive
  # ins = 2 : both p1 and p2 are insertive

  al <- dat$temp$al
# browser()
  # rectal GC infection
  # requires urethral GC in infected partner, infected insertive, and no rectalGC in sus partner
  p1Inf.rgc <- which(urethralGC[al[, "p1"]] == 1 & urethralGC[al[, "p2"]] == 0 &
                     rectalGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf.rgc <- which(urethralGC[al[, "p1"]] == 0 & urethralGC[al[, "p2"]] == 1 &
                     rectalGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs.rgc <- c(p1Inf.rgc, p2Inf.rgc)

  uai.rgc <- al[, "uai"][allActs.rgc]
  patp.rgc <- rep(gc.rt.patp, length(allActs.rgc))
  patp.rgc[uai.rgc == 0] <- patp.rgc[uai.rgc == 0] * sti.cond.rr

  trans.rgc <- rbinom(length(allActs.rgc), 1, patp.rgc)

  transAL.rgc <- al[allActs.rgc[trans.rgc == 1], ]
  idsInf.rgc <- unique(ifelse(urethralGC[transAL.rgc[, "p1"]] == 1,
                              transAL.rgc[, "p2"], transAL.rgc[, "p1"]))

  rectalGC[idsInf.rgc] <- 1

  # urethral GC infection
  # requires rectal GC in infected partner, infected receptive, and no urthethralGC in sus partner
  p1Inf.ugc <- which(rectalGC[al[, "p1"]] == 1 & rectalGC[al[, "p2"]] == 0 &
                     urethralGC[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf.ugc <- which(rectalGC[al[, "p1"]] == 0 & rectalGC[al[, "p2"]] == 1 &
                     urethralGC[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs.ugc <- c(p1Inf.ugc, p2Inf.ugc)

  uai.ugc <- al[, "uai"][allActs.ugc]
  patp.ugc <- rep(gc.ur.patp, length(allActs.ugc))
  patp.ugc[uai.ugc == 0] <- patp.ugc[uai.ugc == 0] * sti.cond.rr

  trans.ugc <- rbinom(length(allActs.ugc), 1, patp.ugc)

  transAL.ugc <- al[allActs.ugc[trans.ugc == 1], ]
  idsInf.ugc <- unique(ifelse(urethralGC[transAL.ugc[, "p1"]] == 1,
                              transAL.ugc[, "p2"], transAL.ugc[, "p1"]))

  urethralGC[idsInf.ugc] <- 1

  # rectal CT infection
  # requires urethral CT in infected partner, infected insertive, and no rectalCT in sus partner
  p1Inf.rct <- which(urethralCT[al[, "p1"]] == 1 & urethralCT[al[, "p2"]] == 0 &
                     rectalCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(1, 2))
  p2Inf.rct <- which(urethralCT[al[, "p1"]] == 0 & urethralCT[al[, "p2"]] == 1 &
                     rectalCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(0, 2))
  allActs.rct <- c(p1Inf.rct, p2Inf.rct)

  uai.rct <- al[, "uai"][allActs.rct]
  patp.rct <- rep(ct.rt.patp, length(allActs.rct))
  patp.rct[uai.rct == 0] <- patp.rct[uai.rct == 0] * sti.cond.rr

  trans.rct <- rbinom(length(allActs.rct), 1, patp.rct)

  transAL.rct <- al[allActs.rct[trans.rct == 1], ]
  idsInf.rct <- unique(ifelse(urethralCT[transAL.rct[, "p1"]] == 1,
                              transAL.rct[, "p2"], transAL.rct[, "p1"]))

  rectalCT[idsInf.rct] <- 1

  # urethral CT infection
  # requires rectal CT in infected partner, infected receptive, and no urthethralCT in sus partner
  p1Inf.uct <- which(rectalCT[al[, "p1"]] == 1 & rectalCT[al[, "p2"]] == 0 &
                       urethralCT[al[, "p2"]] == 0 & al[, "ins"] %in% c(0, 2))
  p2Inf.uct <- which(rectalCT[al[, "p1"]] == 0 & rectalCT[al[, "p2"]] == 1 &
                     urethralCT[al[, "p1"]] == 0 & al[, "ins"] %in% c(1, 2))
  allActs.uct <- c(p1Inf.uct, p2Inf.uct)

  uai.uct <- al[, "uai"][allActs.uct]
  patp.uct <- rep(ct.ur.patp, length(allActs.uct))
  patp.uct[uai.uct == 0] <- patp.uct[uai.uct == 0] * sti.cond.rr

  trans.uct <- rbinom(length(allActs.uct), 1, patp.uct)

  transAL.uct <- al[allActs.uct[trans.uct == 1], ]
  idsInf.uct <- unique(ifelse(urethralCT[transAL.uct[, "p1"]] == 1,
                              transAL.uct[, "p2"], transAL.uct[, "p1"]))

  urethralCT[idsInf.uct] <- 1


  ## Output

  # attributes
  dat$attr$rectalGC <- rectalGC
  dat$attr$urethralGC <- urethralGC
  dat$attr$rectalCT <- rectalCT
  dat$attr$urethralCT <- urethralCT

  # Summary stats
  dat$epi$incid.rgc <- length(idsInf.rgc)
  dat$epi$incid.ugc <- length(idsInf.ugc)
  dat$epi$incid.rct <- length(idsInf.rct)
  dat$epi$incid.uct <- length(idsInf.uct)


  return(dat)
}

dx.sti <- function(dat, at) {

  return(dat)
}

tx.sti <- function(dat, at) {

  return(dat)
}
