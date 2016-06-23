
trans_sti <- function(dat, at) {

  # Variables -----------------------------------------------------------

  # Attributes
  vl <- dat$attr$vl
  stage <- dat$attr$stage
  ccr5 <- dat$attr$ccr5
  circ <- dat$attr$circ
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  tx.status <- dat$attr$tx.status
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  # Parameters
  URAI.prob <- dat$param$URAI.prob
  UIAI.prob <- dat$param$UIAI.prob
  acute.rr <- dat$param$acute.rr
  condom.rr <- dat$param$condom.rr
  circ.rr <- dat$param$circ.rr
  ccr5.heteroz.rr <- dat$param$ccr5.heteroz.rr
  prep.hr <- dat$param$prep.class.hr

  # Data
  al <- dat$temp$al
  dal <- al[which(status[al[, 1]] == 1 & status[al[, 2]] == 0), ]
  dal <- dal[sample(1:nrow(dal)), ]
  ncols <- dim(dal)[2]

  if (nrow(dal) == 0) {
    return(dat)
  }


  ## Reorder by role: ins on the left, rec on the right, flippers represented twice
  disc.ip <- dal[dal[, "ins"] %in% 1:2, ]
  disc.rp <- dal[dal[, "ins"] %in% c(0, 2), c(2:1, 3:ncols)]
  colnames(disc.ip)[1:2] <- colnames(disc.rp)[1:2] <- c("ins", "rec")


  # PATP: Insertive Man Infected (Col 1) --------------------------------

  # Attributes of infected
  ip.vl <- vl[disc.ip[, 1]]
  ip.stage <- stage[disc.ip[, 1]]

  # Attributes of susceptible
  ip.ccr5 <- ccr5[disc.ip[, 2]]
  ip.prep <- prepStat[disc.ip[, 2]]
  ip.prepcl <- prepClass[disc.ip[, 2]]
  ip.rGC <- rGC[disc.ip[, 2]]
  ip.rCT <- rCT[disc.ip[, 2]]

  # Base TP from VL
  trans.ip.prob <- URAI.prob * 2.45^(ip.vl - 4.5)

  # Condom use
  trans.ip.prob[disc.ip[, "uai"] == 0] <- trans.ip.prob[disc.ip[, "uai"] == 0] * condom.rr

  # CCR5
  trans.ip.prob[ip.ccr5 == "DD"] <- trans.ip.prob[ip.ccr5 == "DD"] * 0
  trans.ip.prob[ip.ccr5 == "DW"] <- trans.ip.prob[ip.ccr5 == "DW"] * ccr5.heteroz.rr

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    trans.ip.prob[temp.ids] <- trans.ip.prob[temp.ids] * prep.hr[i]
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c("AR", "AF"))
  trans.ip.prob[isAcute] <- trans.ip.prob[isAcute] * acute.rr

  ## Multiplier for STI
  is.rGC <- which(ip.rGC == 1)
  trans.ip.prob[is.rGC] <- trans.ip.prob[is.rGC] * dat$param$hiv.rgc.rr

  is.rCT <- which(ip.rCT == 1)
  trans.ip.prob[is.rCT] <- trans.ip.prob[is.rCT] * dat$param$hiv.rct.rr


  # PATP: Receptive Man Infected (Col 2) --------------------------------

  # Attributes of infected
  rp.vl <- vl[disc.rp[, 2]]
  rp.stage <- stage[disc.rp[, 2]]

  # Attributes of susceptible
  rp.circ <- circ[disc.rp[, 1]]
  rp.ccr5 <- ccr5[disc.rp[, 1]]
  rp.prep <- prepStat[disc.rp[, 1]]
  rp.prepcl <- prepClass[disc.rp[, 1]]
  rp.uGC <- uGC[disc.rp[, 1]]
  rp.uCT <- uCT[disc.rp[, 1]]

  # Base TP from VL
  trans.rp.prob <- UIAI.prob * 2.45^(rp.vl - 4.5)

  # Circumcision
  trans.rp.prob[rp.circ == 1] <- trans.rp.prob[rp.circ == 1] * circ.rr

  # Condom use
  trans.rp.prob[disc.rp[, "uai"] == 0] <- trans.rp.prob[disc.rp[, "uai"] == 0] * condom.rr

  # CCR5
  trans.rp.prob[rp.ccr5 == "DD"] <- trans.rp.prob[rp.ccr5 == "DD"] * 0
  trans.rp.prob[rp.ccr5 == "DW"] <- trans.rp.prob[rp.ccr5 == "DW"] * ccr5.heteroz.rr

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    trans.rp.prob[temp.ids] <- trans.rp.prob[temp.ids] * prep.hr[i]
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c("AR", "AF"))
  trans.rp.prob[isAcute] <- trans.rp.prob[isAcute] * acute.rr

  ## Multiplier for STI
  is.uGC <- which(rp.uGC == 1)
  trans.rp.prob[is.uGC] <- trans.rp.prob[is.uGC] * dat$param$hiv.ugc.rr

  is.uCT <- which(rp.uCT == 1)
  trans.rp.prob[is.uCT] <- trans.rp.prob[is.uCT] * dat$param$hiv.uct.rr

  ## Bound range of PATP
  trans.ip.prob <- pmin(trans.ip.prob, 1)
  trans.rp.prob <- pmin(trans.rp.prob, 1)

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(trans.ip.prob), 1, trans.ip.prob)
  trans.rp <- rbinom(length(trans.rp.prob), 1, trans.rp.prob)


  ## Output

  # Update attributes

  infected <- infector <- inf.type <- NULL
  if (sum(trans.ip, trans.rp) > 0) {

    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])
    infector <- c(disc.ip[trans.ip == 1, 1],
                  disc.rp[trans.rp == 1, 2])
    inf.role <- c(rep(0, sum(trans.ip)), rep(1, sum(trans.rp)))
    inf.type <- c(disc.ip[trans.ip == 1, "ptype"],
                  disc.rp[trans.rp == 1, "ptype"])

    inf.stage <- stage[infector]
    inf.diag <- diag.status[infector]
    inf.tx <- tx.status[infector]

    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0

    dat$attr$infector[infected] <- infector
    dat$attr$inf.role[infected] <- inf.role
    dat$attr$inf.type[infected] <- inf.type
    dat$attr$inf.diag[infected] <- inf.diag
    dat$attr$inf.tx[infected] <- inf.tx
    dat$attr$inf.stage[infected] <- inf.stage

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
  }

  # Summary Output
  dat$epi$incid[at] <- length(infected)

  if (at >= dat$param$prep.start) {
    dat$epi$mean.trans[at] <- mean(c(trans.ip.prob, trans.rp.prob))
    dat$epi$mean.trans.prep[at] <- mean(c(trans.ip.prob[which(ip.prep == 1)],
                                          trans.rp.prob[which(rp.prep == 1)]))
    dat$epi$mean.trans.nprep[at] <- mean(c(trans.ip.prob[which(ip.prep == 0)],
                                           trans.rp.prob[which(rp.prep == 0)]))
  }

  return(dat)
}
