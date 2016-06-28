
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
  hiv.ugc.rr <- dat$param$hiv.ugc.rr
  hiv.uct.rr <- dat$param$hiv.uct.rr
  hiv.rgc.rr <- dat$param$hiv.rgc.rr
  hiv.rct.rr <- dat$param$hiv.rct.rr
  
  
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
  ip.tprob <- URAI.prob * 2.45^(ip.vl - 4.5)

  # Transform to log odds
  ip.tlo <- log(ip.tprob/(1-ip.tprob))

  # Condom use
  not.UAI <- which(disc.ip[, "uai"] == 0)
  ip.tlo[not.UAI] <- ip.tlo[not.UAI] + log(condom.rr/(1-condom.rr))

  # CCR5
  ip.tlo[ip.ccr5 == "DD"] <- ip.tlo[ip.ccr5 == "DD"] + -Inf
  ip.tlo[ip.ccr5 == "DW"] <- ip.tlo[ip.ccr5 == "DW"] +
                             log(ccr5.heteroz.rr/(1-ccr5.heteroz.rr))

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(ip.prep == 1 & ip.prepcl == i-1)
    ip.tlo[temp.ids] <- ip.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(ip.stage %in% c("AR", "AF"))
  ip.tlo[isAcute] <- ip.tlo[isAcute] + log(acute.rr)

  ## Multiplier for STI
  is.rGC <- which(ip.rGC == 1)
  ip.tlo[is.rGC] <- ip.tlo[is.rGC] + log(hiv.rgc.rr)

  is.rCT <- which(ip.rCT == 1)
  ip.tlo[is.rCT] <- ip.tlo[is.rCT] + log(hiv.rct.rr)

  ip.tprob <- exp(ip.tlo)/(1+exp(ip.tlo))
  stopifnot(ip.tprob >= 0, ip.tprob <= 1)


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
  rp.tprob <- UIAI.prob * 2.45^(rp.vl - 4.5)

  # Transform to log odds
  rp.tlo <- log(rp.tprob/(1-rp.tprob))

  # Circumcision
  rp.tlo[rp.circ == 1] <- rp.tlo[rp.circ == 1] + log(circ.rr/(1-circ.rr))

  # Condom use
  not.UAI <- which(disc.rp[, "uai"] == 0)
  rp.tlo[not.UAI] <- rp.tlo[not.UAI] + log(condom.rr/(1-condom.rr))

  # CCR5
  rp.tlo[rp.ccr5 == "DD"] <- rp.tlo[rp.ccr5 == "DD"] + -Inf
  rp.tlo[rp.ccr5 == "DW"] <- rp.tlo[rp.ccr5 == "DW"] +
                             log(ccr5.heteroz.rr/(1-ccr5.heteroz.rr))

  # PrEP, cycle through 4 adherence classes
  for (i in 1:4) {
    temp.ids <- which(rp.prep == 1 & rp.prepcl == i-1)
    rp.tlo[temp.ids] <- rp.tlo[temp.ids] + log(prep.hr[i])
  }

  # Acute-stage multipliers
  isAcute <- which(rp.stage %in% c("AR", "AF"))
  rp.tlo[isAcute] <- rp.tlo[isAcute] + log(acute.rr)

  ## Multiplier for STI
  is.uGC <- which(rp.uGC == 1)
  rp.tlo[is.uGC] <- rp.tlo[is.uGC] + log(hiv.ugc.rr)

  is.uCT <- which(rp.uCT == 1)
  rp.tlo[is.uCT] <- rp.tlo[is.uCT] + log(hiv.uct.rr)

  # Retransformation to probability
  rp.tprob <- exp(rp.tlo)/(1+exp(rp.tlo))
  stopifnot(rp.tprob >= 0, rp.tprob <= 1)


  # Transmission --------------------------------------------------------

  ## Bernoulli transmission events
  trans.ip <- rbinom(length(ip.tprob), 1, ip.tprob)
  trans.rp <- rbinom(length(rp.tprob), 1, rp.tprob)


  # Output --------------------------------------------------------------


  # Update attributes

  infected <- inf.type <- NULL
  if (sum(trans.ip, trans.rp) > 0) {

    infected <- c(disc.ip[trans.ip == 1, 2],
                  disc.rp[trans.rp == 1, 1])
    inf.role <- c(rep(0, sum(trans.ip)), rep(1, sum(trans.rp)))
    inf.type <- c(disc.ip[trans.ip == 1, "ptype"],
                  disc.rp[trans.rp == 1, "ptype"])

    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0

    dat$attr$inf.role[infected] <- inf.role
    dat$attr$inf.type[infected] <- inf.type

    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
  }

  # Summary Output
  dat$epi$incid[at] <- length(infected)

  return(dat)
}
