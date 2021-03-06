sti_test_msm <- function(dat, at) {

  if (at < dat$param$stitest.start) {
    return(dat)
  }


  ## Variables

  # Attributes
  race <- dat$attr$race
  diag.status.syph <- dat$attr$diag.status.syph
  diag.status.gc <- dat$attr$diag.status.gc
  diag.status.ct <- dat$attr$diag.status.ct
  syphilis <- dat$attr$syphilis
  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  last.neg.test.rgc <- dat$attr$last.neg.test.rgc
  last.neg.test.ugc <- dat$attr$last.neg.test.ugc
  last.neg.test.rct <- dat$attr$last.neg.test.rct
  last.neg.test.uct <- dat$attr$last.neg.test.uct
  last.neg.test.syph <- dat$attr$last.neg.test.syph
  lastdiag.time.gc <- dat$attr$lastdiag.time.gc
  lastdiag.time.ct <- dat$attr$lastdiag.time.ct
  lastdiag.time.syph <- dat$attr$lastdiag.time.syph

  role.class <- dat$attr$role.class
  stage.syph <- dat$attr$stage.syph

  tt.traj.ct <- dat$attr$tt.traj.ct
  tt.traj.gc <- dat$attr$tt.traj.gc
  tt.traj.syph <- dat$attr$tt.traj.syph

  prepStat <- dat$attr$prepStat

  # Parameters
  stianntest.coverage <- dat$param$stianntest.coverage
  stianntest.cov.rate <- dat$param$stianntest.cov.rate
  stihighrisktest.coverage <- dat$param$stihighrisktest.coverage
  stihighrisktest.cov.rate <- dat$param$stihighrisktest.cov.rate
  testing.pattern.sti <- dat$param$testing.pattern.sti
  stitest.active.int <- dat$param$stitest.active.int
  sti.highrisktest.int <- dat$param$sti.highrisktest.int
  tst.rect.sti.rr <- dat$param$tst.rect.sti.rr
  stitest.elig.model <- dat$param$stitest.elig.model

  # Eligibility and trajectory
  # Base eligibility
  idsEligTest <- which(race %in% c("B", "W"))

  # Annual indications- sexually active in last year
  stitestind1 <- dat$attr$stitest.ind.active

  # Annual - testing trajectory update
  activeindwindow <- stitest.active.int * 7
  idsactive <- intersect(which(at - stitestind1 <= activeindwindow), idsEligTest)

  # High-risk - testing trajectory update
  hrindwindow <- sti.highrisktest.int * 7

  #STI testing eligibility scenarios

  # Set to null values if

  if (stitest.elig.model == "all") {
    c1 <- dat$attr$stitest.ind.sti
    c2 <- dat$attr$stitest.ind.recentpartners
    c3 <- dat$attr$stitest.ind.newpartners
    c4 <- dat$attr$stitest.ind.concurrpartner
    c5 <- dat$attr$stitest.ind.partnersti
    c6 <- dat$attr$stitest.ind.uai.nmain
    c7 <- dat$attr$stitest.ind.uai.any

    idshighrisk <- which((at - c1 <= hrindwindow) |
                           (at - c2 <= hrindwindow) |
                           (at - c3 <= hrindwindow) |
                           (at - c4 <= hrindwindow) |
                           (at - c5 <= hrindwindow) |
                           (at - c6 <= hrindwindow) |
                           (at - c7 <= hrindwindow))

    idsnottestelig <- which(tt.traj.syph == 2 & (
      (at - c1 > hrindwindow) & (at - c2 > hrindwindow) &
        (at - c3 > hrindwindow) & (at - c4 > hrindwindow) &
        (at - c5 > hrindwindow) & (at - c6 > hrindwindow) &
        (at - c7 > hrindwindow)))

  } else if (stitest.elig.model != "all") {

    if (substr(stitest.elig.model, 1, 3) == "cdc") {

      if (stitest.elig.model == "cdc1") {
        c1 <- dat$attr$stitest.ind.recentpartners
        c2 <- dat$attr$stitest.ind.newpartners

        idshighrisk <- which((at - c1) <= hrindwindow |
                               (at - c2) <= hrindwindow)

        idsnottestelig <- which(tt.traj.syph == 2 & (
          (at - c1 > hrindwindow) &
            (at - c2 > hrindwindow)))

      } else if (stitest.elig.model == "cdc2") {
        c1 <- dat$attr$stitest.ind.sti
        c2 <- dat$attr$stitest.ind.recentpartners
        c3 <- dat$attr$stitest.ind.newpartners

        idshighrisk <- which((at - c1) <= hrindwindow |
                               (at - c2) <= hrindwindow |
                               (at - c3) <= hrindwindow)

        idsnottestelig <- which(tt.traj.syph == 2 & (
          (at - c1 > hrindwindow) &
            (at - c2 > hrindwindow) &
            (at - c3 > hrindwindow)))

      } else if (stitest.elig.model == "cdc3") {
        c1 <- dat$attr$stitest.ind.uai.nmain
        c2 <- dat$attr$stitest.ind.uai.any

        idshighrisk <- which((at - c1) <= hrindwindow |
                               (at - c2) <= hrindwindow)

        idsnottestelig <- which(tt.traj.syph == 2 & (
          (at - c1 > hrindwindow) &
            (at - c2 > hrindwindow)))

      } else if (stitest.elig.model == "cdc4") {
        c1 <- dat$attr$stitest.ind.concurrpartner
        c2 <- dat$attr$stitest.ind.partnersti

        idshighrisk <- which((at - c1) <= hrindwindow |
                               (at - c2) <= hrindwindow)

        idsnottestelig <- which(tt.traj.syph == 2 & (
          (at - c1 > hrindwindow) &
            (at - c2 > hrindwindow)))

      }
    } else if (stitest.elig.model == "sti") {

      c1 <- dat$attr$stitest.ind.sti

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))

    }  else if (stitest.elig.model == "recentpartners") {

      c1 <- dat$attr$stitest.ind.recentpartners

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))

    }  else if (stitest.elig.model == "newpartners") {

      c1 <- dat$attr$stitest.ind.newpartners

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))
    } else if (stitest.elig.model == "concurrpartner") {

      c1 <- dat$attr$stitest.ind.concurrpartner

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))

    } else if (stitest.elig.model == "partnersti") {

      c1 <- dat$attr$stitest.ind.partnersti

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))

    } else if (stitest.elig.model == "uai.nmain") {

      c1 <- dat$attr$stitest.ind.uai.nmain

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))

    } else if (stitest.elig.model == "uai.any") {

      c1 <- dat$attr$stitest.ind.uai.any

      idshighrisk <- which((at - c1) <= hrindwindow)

      idsnottestelig <- which(tt.traj.syph == 2 &
                                (at - c1 > hrindwindow))

    } else if (stitest.elig.model == "none") {

      idshighrisk <- NULL
      idsnottestelig <- which(race %in% c("B", "W"))

    }
  }

  ## Stoppage (tt.traj.gc/.ct/.syph <- NA------------------------------------
  # Reduce testing trajectory to NA if no longer indicated for more frequent
  # high-risk testing
  dat$attr$stihighrisktestLastElig[idsnottestelig] <- at
  tt.traj.syph[idsnottestelig] <- tt.traj.gc[idsnottestelig] <- tt.traj.ct[idsnottestelig] <- NA

  # Remove testing trajectory if no longer indicated for annual testing
  # (idsannual includes high-risk)
  idsnottestelig <- which(tt.traj.syph == 1 &
                            (at - stitestind1 >= activeindwindow))
  dat$attr$stianntestLastElig[idsnottestelig] <- at
  tt.traj.syph[idsnottestelig] <- tt.traj.gc[idsnottestelig] <-
    tt.traj.ct[idsnottestelig] <- NA


  ## Initiation -------------------------------------------------------------

  ### Testing coverage for high risk
  stihighrisktestCov <- sum(tt.traj.ct == 2, na.rm = TRUE) / length(idshighrisk)
  stihighrisktestCov <- ifelse(is.nan(stihighrisktestCov), 0, stihighrisktestCov)

  idsEligSt <- idshighrisk
  nEligSt <- length(idshighrisk)

  nStart <- max(0, min(nEligSt, round((stihighrisktest.coverage - stihighrisktestCov) *
                                        length(idshighrisk))))
  idsStart <- NULL
  if (nStart > 0) {
    if (stihighrisktest.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, stihighrisktest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory
  if (length(idsStart) > 0) {
    tt.traj.syph[idsStart] <- tt.traj.gc[idsStart] <-tt.traj.ct[idsStart] <- 2
  }

  ### Testing coverage for annual
  # Make this only the sum of where tt.traj.ct == 1? - would need to change
  # denominator to be those who are in not in intersect of idsactive and
  # idshighrisk
  stianntestCov <- sum(tt.traj.ct == 1, na.rm = TRUE) / length(setdiff(idsactive, idshighrisk))
  stianntestCov <- ifelse(is.nan(stianntestCov), 0, stianntestCov)

  idsEligSt <- setdiff(idsactive, idshighrisk)

  nEligSt <- length(setdiff(idsactive, idshighrisk))

  nStart <- max(0, min(nEligSt, round((stianntest.coverage - stianntestCov) *
                                        length(setdiff(idsactive, idshighrisk)))))
  idsStart <- NULL
  if (nStart > 0) {
    if (stianntest.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, stianntest.cov.rate) == 1]
    }
  }

  ## Update testing trajectory
  if (length(idsStart) > 0) {
    tt.traj.syph[idsStart] <- tt.traj.gc[idsStart] <- tt.traj.ct[idsStart] <- 1
  }

  ## Testing
  tsincelntst.syph <- at - dat$attr$last.neg.test.syph
  tsincelntst.syph[is.na(tsincelntst.syph)] <- at - dat$attr$arrival.time[is.na(tsincelntst.syph)]

  tsincelntst.rgc <- at - dat$attr$last.neg.test.rgc
  tsincelntst.ugc <- at - dat$attr$last.neg.test.ugc
  tsincelntst.rgc[is.na(tsincelntst.rgc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rgc)]
  tsincelntst.ugc[is.na(tsincelntst.ugc)] <- at - dat$attr$arrival.time[is.na(tsincelntst.ugc)]
  tsincelntst.gc <- min(tsincelntst.rgc, tsincelntst.ugc)

  tsincelntst.rct <- at - dat$attr$last.neg.test.rct
  tsincelntst.uct <- at - dat$attr$last.neg.test.uct
  tsincelntst.rct[is.na(tsincelntst.rct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.rct)]
  tsincelntst.uct[is.na(tsincelntst.uct)] <- at - dat$attr$arrival.time[is.na(tsincelntst.uct)]
  tsincelntst.ct <- min(tsincelntst.rct, tsincelntst.uct)

  # Testing Rates by serostatus/race?
  # All MSM with HIV infection entering care should be screened for GC and CT
  # ct appropriate anatomic sites of exposure, as well as for syphilis
  # For sexually active individuals, screen at first HIV evaluation,
  # and at least annually thereafter

  # More frequent STD screening (i.e., for syphilis, gonorrhea, and chlamydia)
  # at 3–6-month intervals is indicated for MSM, including those with HIV
  # infection if risk behaviors persist or if they or their sexual partners
  # have multiple partners.

  # Mostly asymptomatic testing handled here - symptomatic testing is equated
  # to probability of symptomatic treatment

  ## Process for syphilis
  if (testing.pattern.sti == "memoryless") {
    elig.syph.ann <- which(tt.traj.syph == 1 &
                             (diag.status.syph == 0 | is.na(diag.status.syph)) &
                             prepStat == 0)
    rates.syph <- rep(1/stitest.active.int, length(elig.syph.ann))
    tst.syph.nprep.ann <- elig.syph.ann[rbinom(length(elig.syph.ann), 1, rates.syph) == 1]

    elig.syph.highrisk <- which(tt.traj.syph == 2 &
                                  (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                  prepStat == 0)
    rates.syph <- rep(1/sti.highrisktest.int, length(elig.syph.highrisk))
    tst.syph.nprep.highrisk <- elig.syph.highrisk[rbinom(length(elig.syph.highrisk), 1, rates.syph) == 1]
    tst.syph.nprep <- c(tst.syph.nprep.ann, tst.syph.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.syph.annual.interval <- which(tt.traj.syph == 1 &
                                        (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                        tsincelntst.syph >= 2*(stitest.active.int) &
                                        prepStat == 0)
    tst.syph.highrisk.interval <- which(tt.traj.syph == 2 &
                                          (diag.status.syph == 0 | is.na(diag.status.syph)) &
                                          tsincelntst.syph >= 2*(sti.highrisktest.int) &
                                          prepStat == 0)
    tst.syph.nprep <- c(tst.syph.annual.interval, tst.syph.highrisk.interval)
  }

  ## Process for GC
  if (testing.pattern.sti == "memoryless") {
    elig.gc.ann <- which(tt.traj.gc == 1 &
                           (diag.status.gc == 0 | is.na(diag.status.gc)) &
                           prepStat == 0)
    rates.gc <- rep(1/stitest.active.int, length(elig.gc.ann))
    tst.gc.nprep.ann <- elig.gc.ann[rbinom(length(elig.gc.ann), 1, rates.gc) == 1]

    elig.gc.highrisk <- which(tt.traj.gc == 2 &
                                (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                prepStat == 0)
    rates.gc <- rep(1/sti.highrisktest.int, length(elig.gc.highrisk))
    tst.gc.nprep.highrisk <- elig.gc.highrisk[rbinom(length(elig.gc.highrisk), 1, rates.gc) == 1]
    tst.gc.nprep <- c(tst.gc.nprep.ann, tst.gc.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.gc.annual.interval <- which(tt.traj.gc == 1 &
                                      (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                      tsincelntst.gc >= 2*(stitest.active.int) &
                                      prepStat == 0)

    tst.gc.highrisk.interval <- which(tt.traj.gc == 2 &
                                        (diag.status.gc == 0 | is.na(diag.status.gc)) &
                                        tsincelntst.gc >= 2*(sti.highrisktest.int) &
                                        prepStat == 0)

    tst.gc.nprep <- c(tst.gc.annual.interval, tst.gc.highrisk.interval)
  }

  ## Process for CT
  if (testing.pattern.sti == "memoryless") {
    elig.ct.ann <- which(tt.traj.ct == 1 &
                           (diag.status.ct == 0 | is.na(diag.status.ct)) &
                           prepStat == 0)
    rates.ct <- rep(1/stitest.active.int, length(elig.ct.ann))
    tst.ct.nprep.ann <- elig.ct.ann[rbinom(length(elig.ct.ann), 1, rates.ct) == 1]

    elig.ct.highrisk <- which(tt.traj.ct == 2 &
                                (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                prepStat == 0)
    rates.ct <- rep(1/sti.highrisktest.int, length(elig.ct.highrisk))
    tst.ct.nprep.highrisk <- elig.ct.highrisk[rbinom(length(elig.ct.highrisk), 1, rates.ct) == 1]

    tst.ct.nprep <- c(tst.ct.nprep.ann, tst.ct.nprep.highrisk)
  }

  if (testing.pattern.sti == "interval" ) {
    tst.ct.annual.interval <- which(tt.traj.ct == 1 &
                                      (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                      tsincelntst.ct >= 2*(stitest.active.int) &
                                      prepStat == 0)

    tst.ct.highrisk.interval <- which(tt.traj.ct == 2 &
                                        (diag.status.ct == 0 | is.na(diag.status.ct)) &
                                        tsincelntst.ct >= 2*(sti.highrisktest.int) &
                                        prepStat == 0)

    tst.ct.nprep <- c(tst.ct.annual.interval, tst.ct.highrisk.interval)
  }

  # Syphilis non-PrEP testing
  tst.syph.pos <- tst.syph.nprep[syphilis[tst.syph.nprep] == 1 &
                                   stage.syph[tst.syph.nprep] %in% c(2, 3, 4, 5, 6, 7)]
  tst.syph.neg <- setdiff(tst.syph.nprep, tst.syph.pos)

  # GC non-PrEP testing
  tst.rgc <- tst.gc.nprep[role.class[tst.gc.nprep] %in% c("R", "V")]
  tst.rgc <- sample(tst.rgc, tst.rect.sti.rr * length(tst.rgc))
  tst.ugc <- tst.gc.nprep[role.class[tst.gc.nprep] %in% c("I", "V")]
  tst.rgc.pos <- tst.rgc[rGC[tst.rgc] == 1]
  tst.ugc.pos <- tst.ugc[uGC[tst.ugc] == 1]
  tst.rgc.neg <- setdiff(tst.rgc, tst.rgc.pos)
  tst.ugc.neg <- setdiff(tst.ugc, tst.ugc.pos)
  tst.gc.pos <- unique(c(tst.rgc.pos, tst.ugc.pos))

  # CT non-PrEP testing
  tst.rct <- tst.ct.nprep[role.class[tst.ct.nprep] %in% c("R", "V")]
  tst.rct <- sample(tst.rct, tst.rect.sti.rr * length(tst.rct))
  tst.uct <- tst.ct.nprep[role.class[tst.ct.nprep] %in% c("I", "V")]
  tst.rct.pos <- tst.rct[rCT[tst.rct] == 1]
  tst.uct.pos <- tst.uct[uCT[tst.uct] == 1]
  tst.rct.neg <- setdiff(tst.rct, tst.rct.pos)
  tst.uct.neg <- setdiff(tst.uct, tst.uct.pos)
  tst.ct.pos <- unique(c(tst.rct.pos, tst.uct.pos))

  # Syphilis Attributes
  last.neg.test.syph[tst.syph.neg] <- at
  last.neg.test.syph[tst.syph.pos] <- NA
  diag.status.syph[tst.syph.pos] <- 1
  lastdiag.time.syph[tst.syph.pos] <- at

  # GC Attributes
  last.neg.test.rgc[tst.rgc.neg] <- at
  last.neg.test.ugc[tst.ugc.neg] <- at
  last.neg.test.rgc[tst.rgc.pos] <- NA
  last.neg.test.ugc[tst.ugc.pos] <- NA
  diag.status.gc[tst.gc.pos] <- 1
  lastdiag.time.gc[tst.gc.pos] <- at

  # CT Attributes
  last.neg.test.rct[tst.rct.neg] <- at
  last.neg.test.uct[tst.uct.neg] <- at
  last.neg.test.rct[tst.rct.pos] <- NA
  last.neg.test.uct[tst.uct.pos] <- NA
  diag.status.ct[tst.ct.pos] <- 1
  lastdiag.time.ct[tst.ct.pos] <- at

  if (is.null(dat$epi$num.asympt.tx)) {
    dat$epi$rGCasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$uGCasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$GCasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$rCTasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$uCTasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$CTasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$syphasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$totalstiasympttests <- rep(0, length(dat$control$nsteps))
    dat$epi$rGCasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$uGCasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$GCasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$rCTasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$uCTasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$CTasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$syphasympttests.pos <- rep(0, length(dat$control$nsteps))
    dat$epi$totalstiasympttests.pos <- rep(0, length(dat$control$nsteps))
  }

  # Number of tests for asymptomatic
  dat$epi$rGCasympttests[at] <- length(tst.rgc)
  dat$epi$uGCasympttests[at] <- length(tst.ugc)
  dat$epi$GCasympttests[at] <- length(c(tst.rgc, tst.ugc))

  dat$epi$rGCasympttests.pos[at] <- length(tst.rgc.pos)
  dat$epi$uGCasympttests.pos[at] <- length(tst.ugc.pos)
  dat$epi$GCasympttests.pos[at] <- length(c(tst.rgc.pos, tst.ugc.pos))

  dat$epi$rCTasympttests[at] <- length(tst.rct)
  dat$epi$uCTasympttests[at] <- length(tst.uct)
  dat$epi$CTasympttests[at] <- length(c(tst.rct, tst.uct))

  dat$epi$rCTasympttests.pos[at] <- length(tst.rct.pos)
  dat$epi$uCTasympttests.pos[at] <- length(tst.uct.pos)
  dat$epi$CTasympttests.pos[at] <- length(c(tst.rct.pos, tst.uct.pos))

  dat$epi$syphasympttests[at] <- length(c(tst.syph.nprep))
  dat$epi$syphasympttests.pos[at] <- length(c(tst.syph.pos))

  dat$epi$totalstiasympttests[at] <- length(c(tst.rct, tst.uct, tst.rgc,
                                              tst.ugc, tst.syph.nprep))
  dat$epi$totalstiasympttests.pos[at] <- length(c(tst.rgc.pos, tst.ugc.pos,
                                                  tst.rct.pos, tst.uct.pos,
                                                  tst.syph.pos))


  ## Output -----------------------------------------------------------------

  # Attributes

  # Syphilis Attributes
  dat$attr$last.neg.test.syph <- last.neg.test.syph
  dat$attr$diag.status.syph <- diag.status.syph
  dat$attr$lastdiag.time.syph <- lastdiag.time.syph
  dat$attr$tt.traj.syph <- tt.traj.syph

  # GC Attributes
  dat$attr$last.neg.test.rgc <- last.neg.test.rgc
  dat$attr$last.neg.test.ugc <- last.neg.test.ugc
  dat$attr$diag.status.gc <- diag.status.gc
  dat$attr$lastdiag.time.gc <- lastdiag.time.gc
  dat$attr$tt.traj.gc <- tt.traj.gc

  # CT Attributes
  dat$attr$last.neg.test.rct <- last.neg.test.rct
  dat$attr$last.neg.test.uct <- last.neg.test.uct
  dat$attr$diag.status.ct <- diag.status.ct
  dat$attr$lastdiag.time.ct <- lastdiag.time.ct
  dat$attr$tt.traj.ct <- tt.traj.ct


  return(dat)
}


#' @title HIV Diagnosis Module
#'
#' @description Module function for simulating HIV diagnosis after infection,
#'              currently based on diagnosis at treatment initiation.
#'
#' @inheritParams aging_het
#'
#' @keywords module het
#'
#' @export
#'
dx_het <- function(dat, at) {

  # Variables
  status <- dat$attr$status
  txCD4min <- dat$attr$txCD4min
  cd4Count <- dat$attr$cd4Count
  dxStat <- dat$attr$dxStat

  # Process
  tested <- which(status == 1 & dxStat == 0 & cd4Count <= txCD4min)


  # Results
  if (length(tested) > 0) {
    dat$attr$dxStat[tested] <- 1
    dat$attr$txStat[tested] <- 0
    dat$attr$dxTime[tested] <- at
  }

  return(dat)
}
