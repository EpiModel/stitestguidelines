
### prep.FUN = prep.mard

prep.sti <- function(dat, at) {

  if (at < dat$param$prep.start) {
    return(dat)
  }

  ## Variables
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk

  prep.coverage <- dat$param$prep.coverage
  prep.cov.method <- dat$param$prep.cov.method
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob


  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0 & lnt == at)

  # Core eligiblity
  mat.c1 <- dat$riskh$uai.mono
  mat.c2 <- dat$riskh$uai.nmain
  mat.c3 <- dat$riskh$ai.sd
  mat.c4 <- dat$riskh$sti

  idsEligStart <- intersect(which(rowSums(mat.c1, na.rm = TRUE) > 0 |
                                  rowSums(mat.c2, na.rm = TRUE) > 0 |
                                  rowSums(mat.c3, na.rm = TRUE) > 0 |
                                  rowSums(mat.c4, na.rm = TRUE) > 0),
                            idsEligStart)

  prepElig[idsEligStart] <- 1


  ## Stoppage ------------------------------------------------------------------

  # No indications
  idsRiskAssess <- which(active == 1 & prepStat == 1 & lnt == at & (at - prepLastRisk) >= 52)
  prepLastRisk[idsRiskAssess] <- at

  idsEligStop <- intersect(which(rowSums(mat.c1, na.rm = TRUE) == 0 &
                                   rowSums(mat.c2, na.rm = TRUE) == 0 &
                                   rowSums(mat.c3, na.rm = TRUE) == 0 &
                                   rowSums(mat.c4, na.rm = TRUE) == 0),
                           idsRiskAssess)
  prepElig[idsEligStop] <- 0

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsEligStop)
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  if (prep.cov.method == "curr") {
    prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov <- sum(prepEver == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(prepElig == 1)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                        sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0) {
    if (prep.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, prep.cov.rate) == 1]
    }
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepEver[idsStart] <- 1
    prepLastRisk[idsStart] <- at

    # PrEP class
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = 0:3, size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass
  dat$attr$prepLastRisk <- prepLastRisk

  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)

  return(dat)
}
