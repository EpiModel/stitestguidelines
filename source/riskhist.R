
### riskhist.FUN = riskhist.mard

riskhist.sti <- function(dat, at) {

  if (at < dat$param$riskh.start) {
    return(dat)
  }

  ## Attributes
  uid <- dat$attr$uid

  ## Parameters
  pri <- ceiling(dat$param$prep.risk.int)

  ## Edgelist, adds uai summation per partnership from act list
  al <- dat$temp$al
  uai <- as.numeric(by(al[, "uai"], al[, "pid"], sum))
  el <- as.data.frame(cbind(dat$temp$el, uai))

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  ## Truncate riskh matrices
  for (i in 1:length(dat$riskh)) {
    nc <- ncol(dat$riskh[[i]])
    if (pri < ncol(dat$riskh[[i]])) {
      dat$riskh[[i]] <- dat$riskh[[i]][, (nc - pri + 1):nc]
    }
    if (pri > nc) {
      nr <- nrow(dat$riskh[[i]])
      dat$riskh[[i]] <- cbind(matrix(NA, ncol = (pri - nc), nrow = nr),
                              dat$riskh[[i]])
    }
    dat$riskh[[i]] <- dat$riskh[[i]][, -1]
    dat$riskh[[i]] <- cbind(dat$riskh[[i]], rep(NA, nrow(dat$riskh[[i]])))
  }

  ## Degree ##
  n <- attributes(dat$el[[1]])$n
  main.deg <- casl.deg <- inst.deg <- rep(0, n)

  tab.main <- table(dat$el[[1]])
  main.deg[as.numeric(names(tab.main))] <- as.vector(tab.main)

  tab.casl <- table(dat$el[[2]])
  casl.deg[as.numeric(names(tab.casl))] <- as.vector(tab.casl)

  tab.inst <- table(dat$el[[3]])
  inst.deg[as.numeric(names(tab.inst))] <- as.vector(tab.inst)


  ## Preconditions ##

  # Any UAI
  uai.any <- unique(c(el2$p1[el2$uai > 0],
                      el2$p2[el2$uai > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg + casl.deg + inst.deg
  uai.mono1 <- intersect(which(tot.deg == 1), uai.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  dx <- dat$attr$diag.status
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)
  since.test <- at - dat$attr$last.neg.test

  ## Condition 1b: UAI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  uai.mono1.neg <- intersect(uai.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% uai.mono1.neg, 2], el2[el2$p2 %in% uai.mono1.neg, 1])
  not.tested.6mo <- since.test[part.id1] > (180/dat$param$time.unit)
  part.not.tested.6mo <- uai.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$riskh$uai.mono[, pri] <- 0
  dat$riskh$uai.mono[part.not.tested.6mo, pri] <- 1

  ## Condition 2b: UAI in non-main partnerships
  uai.nmain <- unique(c(el2$p1[el2$st1 == 0 & el2$uai > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$uai > 0 & el2$ptype %in% 2:3]))
  dat$riskh$uai.nmain[, pri] <- 0
  dat$riskh$uai.nmain[uai.nmain, pri] <- 1

  ## Condition 3a: AI within known serodiscordant partnerships
  el2.cond3 <- el2[el2$st1 == 1 & el2$ptype %in% 1:2, ]

  # Disclosure
  discl.list <- dat$temp$discl.list
  disclose.cdl <- discl.list[, 1] * 1e7 + discl.list[, 2]
  delt.cdl <- uid[el2.cond3[, 1]] * 1e7 + uid[el2.cond3[, 2]]
  discl <- (delt.cdl %in% disclose.cdl)

  ai.sd <- el2.cond3$p2[discl == TRUE]
  dat$riskh$ai.sd[, pri] <- 0
  dat$riskh$ai.sd[ai.sd, pri] <- 1


  ## Condition 4, any STI diagnosis
  idsDx <- which(dat$attr$rGC.tx == 1 | dat$attr$uGC.tx == 1 |
                 dat$attr$rCT.tx == 1 | dat$attr$uCT.tx == 1)
  dat$riskh$sti[, pri] <- 0
  dat$riskh$sti[idsDx, pri] <- 1

  return(dat)
}

