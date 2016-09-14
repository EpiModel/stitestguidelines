
## P2 functions

gc_pia_vec <- function(sim.base, sim.comp) {

  sim.base <- truncate_sim(sim.base, at = 2600)
  sim.comp <- truncate_sim(sim.comp, at = 2600)

  ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000

  ir.comp.gc <- unname(colMeans(sim.comp$epi$ir100.gc, na.rm = TRUE)) * 1000
  vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)

  vec.pia.gc <- vec.nia.gc/ir.base.gc

  return(vec.pia.gc)
}


epi_stats <- function(sim.base,
                      sim.comp = NULL,
                      at = 520,
                      qnt.low = 0.025,
                      qnt.high = 0.975) {

  # Base scenario -------------------------------------------------------

  sim.base <- truncate_sim(sim.base, at = 2600)

  mn <- as.data.frame(sim.base)
  ql <- as.data.frame(sim.base, out = "qnt", qval = qnt.low)
  qu <- as.data.frame(sim.base, out = "qnt", qval = qnt.high)

  # prevalence
  prev.base <- round(data.frame(mean = mn$i.prev[at],
                                ql = ql$i.prev[at],
                                qu = qu$i.prev[at]), 3)

  # incidence
  haz.base <- round(data.frame(mean = mean(tail(mn$ir100, 52)),
                               ql = mean(tail(ql$ir100, 52)),
                               qu = mean(tail(qu$ir100, 52))), 2)
  ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
  incid.base <- unname(colSums(sim.base$epi$incid))

  haz.base.gc <- round(data.frame(mean = mean(tail(mn$ir100.gc, 52)),
                                  ql = mean(tail(ql$ir100.gc, 52)),
                                  qu = mean(tail(qu$ir100.gc, 52))), 2)
  ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000

  haz.base.ct <- round(data.frame(mean = mean(tail(mn$ir100.ct, 52)),
                                  ql = mean(tail(ql$ir100.ct, 52)),
                                  qu = mean(tail(qu$ir100.ct, 52))), 2)
  ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000


  # Comparison scenario -------------------------------------------------


  if (!is.null(sim.comp)) {

    sim.comp <- truncate_sim(sim.comp, at = 2600)

    mn.comp <- as.data.frame(sim.comp)
    ql.comp <- as.data.frame(sim.comp, out = "qnt", qval = qnt.low)
    qu.comp <- as.data.frame(sim.comp, out = "qnt", qval = qnt.high)

    # # prevalence
    # out.prev <- round(data.frame(mean = mn.comp$i.prev[at],
    #                               ql = ql.comp$i.prev[at],
    #                               qu = qu.comp$i.prev[at]), 3)

    # incidence
    out.haz <- round(data.frame(mean = mean(tail(mn.comp$ir100, 52)),
                                ql = mean(tail(ql.comp$ir100, 52)),
                                qu = mean(tail(qu.comp$ir100, 52))), 2)
    out.haz.gc <- round(data.frame(mean = mean(tail(mn.comp$ir100.gc, 52)),
                                   ql = mean(tail(ql.comp$ir100.gc, 52)),
                                   qu = mean(tail(qu.comp$ir100.gc, 52))), 2)
    out.haz.ct <- round(data.frame(mean = mean(tail(mn.comp$ir100.ct, 52)),
                                   ql = mean(tail(ql.comp$ir100.ct, 52)),
                                   qu = mean(tail(qu.comp$ir100.ct, 52))), 2)

    # HR
    num <- unname(colMeans(tail(sim.comp$epi$ir100, 52)))
    denom <- unname(colMeans(tail(sim.base$epi$ir100, 52)))
    vec.hr <- num/denom
    out.hr <- round(data.frame(mean = mean(vec.hr),
                               ql = quantile(vec.hr, qnt.low, names = FALSE),
                               qu = quantile(vec.hr, qnt.high, names = FALSE)), 3)

    num.gc <- unname(colMeans(tail(sim.comp$epi$ir100.gc, 52)))
    denom.gc <- unname(colMeans(tail(sim.base$epi$ir100.gc, 52)))
    vec.hr.gc <- num.gc/denom.gc
    out.hr.gc <- round(data.frame(mean = mean(vec.hr.gc),
                                  ql = quantile(vec.hr.gc, qnt.low, names = FALSE),
                                  qu = quantile(vec.hr.gc, qnt.high, names = FALSE)), 3)

    num.ct <- unname(colMeans(tail(sim.comp$epi$ir100.ct, 52)))
    denom.ct <- unname(colMeans(tail(sim.base$epi$ir100.ct, 52)))
    vec.hr.ct <- num.ct/denom.ct
    out.hr.ct <- round(data.frame(mean = mean(vec.hr.ct),
                                  ql = quantile(vec.hr.ct, qnt.low, names = FALSE),
                                  qu = quantile(vec.hr.ct, qnt.high, names = FALSE)), 3)

    # NIA
    ir.comp <- unname(colMeans(sim.comp$epi$ir100)) * 1000
    vec.nia <- round(ir.base - ir.comp, 1)
    out.nia <- round(data.frame(mean = mean(vec.nia),
                                ql = quantile(vec.nia, qnt.low, names = FALSE),
                                qu = quantile(vec.nia, qnt.high, names = FALSE)), 0)

    ir.comp.gc <- unname(colMeans(sim.comp$epi$ir100.gc)) * 1000
    vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)

    ir.comp.ct <- unname(colMeans(sim.comp$epi$ir100.ct)) * 1000
    vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)

    # PIA
    vec.pia <- vec.nia/ir.base
    out.pia <- round(data.frame(mean = mean(vec.pia),
                                ql = quantile(vec.pia, qnt.low, names = FALSE),
                                qu = quantile(vec.pia, qnt.high, names = FALSE)), 3)

    vec.pia.gc <- vec.nia.gc/ir.base.gc
    out.pia.gc <- round(data.frame(mean = mean(vec.pia.gc),
                                   ql = quantile(vec.pia.gc, qnt.low, names = FALSE),
                                   qu = quantile(vec.pia.gc, qnt.high, names = FALSE)), 3)

    vec.pia.ct <- vec.nia.ct/ir.base.ct
    out.pia.ct <- round(data.frame(mean = mean(vec.pia.ct),
                                   ql = quantile(vec.pia.ct, qnt.low, names = FALSE),
                                   qu = quantile(vec.pia.ct, qnt.high, names = FALSE)), 3)

    # NNT
    py.on.prep <- unname(colSums(sim.comp$epi$prepCurr))/52
    vec.nnt <- py.on.prep/(incid.base - unname(colSums(sim.comp$epi$incid)))
    out.nnt <- round(data.frame(mean = mean(vec.nnt),
                                ql = quantile(vec.nnt, qnt.low, names = FALSE),
                                qu = quantile(vec.nnt, qnt.high, names = FALSE)), 0)

    # cat("\n\n\Prevalence:")
    # print(t(out.prev))

    cat("\nHIV Incidence:")
    print(t(out.haz))

    cat("\nHIV HR:")
    print(t(out.hr))

    # cat("\nHIV NIA:")
    # print(t(out.nia))

    cat("\nHIV PIA:")
    print(t(out.pia))

    # cat("\nNNT:")
    # print(t(out.nnt))

    cat("\nGC Incidence:")
    print(t(out.haz.gc))

    cat("\nGC HR:")
    print(t(out.hr.gc))

    cat("\nGC PIA:")
    print(t(out.pia.gc))

    cat("\nCT Incidence:")
    print(t(out.haz.ct))

    cat("\nCT HR:")
    print(t(out.hr.ct))

    cat("\nCT PIA:")
    print(t(out.pia.ct))

  } else {

    # cat("\n\nPrevalence:")
    # print(t(prev.base))

    cat("\nHIV Incidence:")
    print(t(haz.base))

    cat("\nGC Incidence:")
    print(t(haz.base.gc))

    cat("\nCT Incidence:")
    print(t(haz.base.ct))

  }

}
