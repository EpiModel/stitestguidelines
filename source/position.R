
position_sti <- function(dat, at) {

  ## Variables
  al <- dat$temp$al
  if (nrow(al) == 0) {
    return(dat)
  }

  # Attributes
  
  role.class <- dat$attr$role.class
  ins.quot <- dat$attr$ins.quot
  race <- dat$attr$race

  # Parameters
  
  vv.iev.BB.prob <- dat$param$vv.iev.BB.prob
  vv.iev.BW.prob <- dat$param$vv.iev.BW.prob
  vv.iev.WW.prob <- dat$param$vv.iev.WW.prob


  ## Process
  p1.role.class <- role.class[al[, 1]]
  p2.role.class <- role.class[al[, 2]]

  ins <- rep(NA, length(p1.role.class))
  ins[which(p1.role.class == "I")] <- 1
  ins[which(p1.role.class == "R")] <- 0
  ins[which(p2.role.class == "I")] <- 0
  ins[which(p2.role.class == "R")] <- 1

  vv <- which(p1.role.class == "V" & p2.role.class == "V")
  vv.race.combo <- paste0(race[al[, 1]][vv], race[al[, 2]][vv])
  vv.race.combo[vv.race.combo == "WB"] <- "BW"
  vv.iev.prob <- (vv.race.combo == "BB") * vv.iev.BB.prob +
                 (vv.race.combo == "BW") * vv.iev.BW.prob +
                 (vv.race.combo == "WW") * vv.iev.WW.prob

  # intra-event versatility
  iev <- rbinom(length(vv), 1, vv.iev.prob)
  ins[vv[iev == 1]] <- 2 # both are insertive, acts will be doubled
  vv.remaining <- vv[iev == 0]

  p1.ins.prob <- ins.quot[al[, 1][vv.remaining]] /
                 (ins.quot[al[, 1][vv.remaining]] + ins.quot[al[, 2][vv.remaining]])
  p1.ins <- rbinom(length(vv.remaining), 1, p1.ins.prob)
  ins[vv.remaining[p1.ins == 1]] <- 1
  ins[vv.remaining[p1.ins == 0]] <- 0


  ## Output
  dat$temp$al <- cbind(al, ins)

  return(dat)
}

