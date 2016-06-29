
initialize_sti <- function(x, param, init, control, s) {

  ## Variables

  # Attributes
  dat$attr$rGC <- rGC
  dat$attr$uGC <- uGC
    
  dat$attr$rCT <- rCT
  dat$attr$uCT <- uCT
  sqrt.age <- dat$attr$sqrt.age
  
  # Parameters
  
    
  
  # Master data list
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control

  dat$attr <- list()
  dat$stats <- list()
  dat$stats$nwstats <- list()
  dat$temp <- list()
  dat$epi <- list()

  ## Network simulation ##
  nw <- list()
  for (i in 1:3) {
    nw[[i]] <- simulate(x[[i]]$fit)
    nw[[i]] <- remove_bad_roles_msm(nw[[i]])
  }

  ## Build initial edgelists
  dat$el <- list()
  dat$p <- list()
  for (i in 1:2) {
    dat$el[[i]] <- as.edgelist(nw[[i]])
    attributes(dat$el[[i]])$vnames <- NULL
    p <- tergmLite::stergm_prep(nw[[i]], x[[i]]$formation, x[[i]]$coef.diss$dissolution,
                                x[[i]]$coef.form, x[[i]]$coef.diss$coef.adj, x[[i]]$constraints)
    p$model.form$formula <- NULL
    p$model.diss$formula <- NULL
    dat$p[[i]] <- p
  }
  dat$el[[3]] <- as.edgelist(nw[[3]])
  attributes(dat$el[[3]])$vnames <- NULL
  p <- tergmLite::ergm_prep(nw[[3]], x[[3]]$formation, x[[3]]$coef.form, x[[3]]$constraints)
  p$model.form$formula <- NULL
  dat$p[[3]] <- p


  # Network parameters
  dat$nwparam <- list()
  for (i in 1:3) {
    dat$nwparam[i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }


  ## Nodal attributes ##

  # Degree terms
  dat$attr$deg.pers <- get.vertex.attribute(x[[1]]$fit$network, "deg.pers")
  dat$attr$deg.main <- get.vertex.attribute(x[[2]]$fit$network, "deg.main")


  # Race
  dat$attr$race <- get.vertex.attribute(nw[[1]], "race")
  num.B <- dat$init$num.B
  num.W <- dat$init$num.W
  num <- num.B + num.W
  ids.B <- which(dat$attr$race == "B")
  ids.W <- which(dat$attr$race == "W")

  dat$attr$active <- rep(1, num)
  dat$attr$uid <- 1:num
  dat$temp$max.uid <- num

  # Age
  dat$attr$sqrt.age <- get.vertex.attribute(nw[[1]], "sqrt.age")
  dat$attr$age <- sqrt.age^2

  # Risk group
  dat$attr$riskg <- get.vertex.attribute(nw[[3]], "riskg")

  # UAI group
  p1 <- dat$param$cond.pers.always.prob
  p2 <- dat$param$cond.inst.always.prob
  rho <- dat$param$cond.always.prob.corr
  uai.always <- bindata::rmvbin(num, c(p1, p2), bincorr = (1 - rho) * diag(2) + rho)
  dat$attr$cond.always.pers <- uai.always[, 1]
  dat$attr$cond.always.inst <- uai.always[, 2]

  # Arrival and departure
  dat$attr$arrival.time <- rep(1, num)

  # Circumcision
  circ <- rep(NA, num)
  circ[ids.B] <- sample(apportion_lr(num.B, 0:1, 1 - param$circ.B.prob))
  circ[ids.W] <- sample(apportion_lr(num.W, 0:1, 1 - param$circ.W.prob))
  dat$attr$circ <- circ

  # PrEP Attributes
  dat$attr$prepClass <- rep(NA, num)
  dat$attr$prepElig <- rep(NA, num)
  dat$attr$prepStat <- rep(0, num)
  dat$attr$prepStartTime <- rep(NA, num)
  dat$attr$prepLastRisk <- rep(NA, num)
  dat$attr$prepLastStiScreen <- rep(NA, num)

  # Risk history lists
  nc <- ceiling(dat$param$prep.risk.int)
  dat$riskh <- list()
  rh.names <- c("uai.mono", "uai.nmain", "ai.sd", "sti")
  for (i in 1:length(rh.names)) {
    dat$riskh[[rh.names[i]]] <- matrix(NA, ncol = nc, nrow = num)
  }


  # One-off AI class
  inst.ai.class <- rep(NA, num)
  ncl <- param$num.inst.ai.classes
  inst.ai.class[ids.B] <- sample(apportion_lr(num.B, 1:ncl, rep(1 / ncl, ncl)))
  inst.ai.class[ids.W] <- sample(apportion_lr(num.W, 1:ncl, rep(1 / ncl, ncl)))
  dat$attr$inst.ai.class <- inst.ai.class

  # Role class
  role.class <- network::get.vertex.attribute(nw[[1]], "role.class")
  dat$attr$role.class <- role.class

  # Ins.quot
  ins.quot <- rep(NA, num)
  ins.quot[role.class == "I"]  <- 1
  ins.quot[role.class == "R"]  <- 0
  ins.quot[role.class == "V"]  <- runif(sum(role.class == "V"))
  dat$attr$ins.quot <- ins.quot

  # HIV-related attributes
  dat <- init_status_msm(dat)

  ## GC/CT status
  idsUreth <- which(role.class %in% c("I", "V"))
  idsRect <- which(role.class %in% c("R", "V"))

  uGC <- rGC <- rep(0, num)
  uCT <- rCT <- rep(0, num)

  # Initialize GC infection at both sites
  idsUGC <- sample(idsUreth, size = round(init$prev.ugc * num), FALSE)
  uGC[idsUGC] <- 1

  idsRGC <- sample(setdiff(idsRect, idsUGC), size = round(init$prev.rgc * num), FALSE)
  rGC[idsRGC] <- 1


  dat$attr$rGC.sympt <- dat$attr$uGC.sympt <- rep(NA, num)
  dat$attr$rGC.sympt[rGC == 1] <- rbinom(sum(rGC == 1), 1, dat$param$rgc.sympt.prob)
  dat$attr$uGC.sympt[uGC == 1] <- rbinom(sum(uGC == 1), 1, dat$param$ugc.sympt.prob)

  dat$attr$rGC.infTime <- dat$attr$uGC.infTime <- rep(NA, length(dat$attr$active))
  dat$attr$rGC.infTime[rGC == 1] <- 1
  dat$attr$uGC.infTime[uGC == 1] <- 1

  dat$attr$rGC.tx <- dat$attr$uGC.tx <- rep(NA, num)
  dat$attr$GC.cease <- rep(NA, num)

  # Initialize CT infection at both sites
  idsUCT <- sample(idsUreth, size = round(init$prev.uct * num), FALSE)
  uCT[idsUCT] <- 1

  idsRCT <- sample(setdiff(idsRect, idsUCT), size = round(init$prev.rct * num), FALSE)
  rCT[idsRCT] <- 1


  dat$attr$rCT.sympt <- dat$attr$uCT.sympt <- rep(NA, num)
  dat$attr$rCT.sympt[rCT == 1] <- rbinom(sum(rCT == 1), 1, dat$param$rct.sympt.prob)
  dat$attr$uCT.sympt[uCT == 1] <- rbinom(sum(uCT == 1), 1, dat$param$uct.sympt.prob)

  dat$attr$rCT.infTime <- dat$attr$uCT.infTime <- rep(NA, num)
  dat$attr$rCT.infTime[dat$attr$rCT == 1] <- 1
  dat$attr$uCT.infTime[dat$attr$uCT == 1] <- 1

  dat$attr$rCT.tx <- dat$attr$uCT.tx <- rep(NA, num)
  dat$attr$CT.cease <- rep(NA, num)

  # CCR5
  dat <- init_ccr5_msm(dat)


  # Network statistics
  dat$stats$nwstats <- list()


  # Prevalence Tracking
  dat$temp$deg.dists <- list()
  dat$temp$discl.list <- matrix(NA, nrow = 0, ncol = 3)
  colnames(dat$temp$discl.list) <- c("pos", "neg", "discl.time")

  dat <- prevalence_sti(dat, at = 1)

  class(dat) <- "dat"
  return(dat)
}
