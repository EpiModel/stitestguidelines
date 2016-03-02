
## STI Modules for stiPrEP Model

## TODO

sti_trans <- function(dat, at) {

  browser()

  ## Parameters
  gc.rt.patp <- dat$param$gc.rt.patp
  gc.ur.patp <- dat$param$gc.ur.patp
  ct.rt.patp <- dat$param$ct.rt.patp
  ct.ur.patp <- dat$param.ct.ur.patp

  sti.cond.rr <- dat$param$sti.cond.rr

  ## Attributes
  rectalGC <- dat$attr$rectalGC
  urethralGC <- dat$attr$urethralGC
  rectalCT <- dat$attr$rectalCT
  urethralCT <- dat$attr$urethralCT

  ## Processes

  # rectal GC infection


  # urethral GC infection


  # rectal CT infection


  # urethral CT infection


  ## Summary statistics



  ## Output


  return(dat)
}

dx.sti <- function(dat, at) {

  return(dat)
}

tx.sti <- function(dat, at) {

  return(dat)
}
