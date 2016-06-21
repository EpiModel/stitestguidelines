
verbose_sti <- function(x, type, s, at) {

  if (type == "startup") {

  }

  if (type == "progress") {
    if (x$control$verbose == TRUE) {
      if (!interactive()) {
        if (!is.null(x$control$ncores) && x$control$ncores > 1) {
          if (x$control$verbose.int > 0 && (at %% x$control$verbose.int == 0)) {
            simno <- x$control$simno
            currsim <- x$control$currsim
            if (file.exists("verb/") == FALSE) {
              dir.create("verb/")
            }
            fn <- paste0("verb/sim", simno, ".s", currsim, ".txt")
            cat("SIMNO ", paste(simno, currsim, sep = "."),
                "\n====================",
                "\nStep: ", at, " (", round(at/x$control$nsteps, 2), ")",
                "\nPop Size: ", x$epi$num[at],
                "\nTot Prev: ", round(x$epi$i.num[at] / x$epi$num[at], 3),
                "\n\n", sep = "", file = fn)
          }
        }
      } else {
        verbose.int <- x$control$verbose.int
        if (verbose.int > 0 && (at %% verbose.int == 0)) {

          nsteps <- x$control$nsteps
          time.unit <- x$param$time.unit
          prev <- round(x$epi$i.prev[at], 3)
          prev.rgc <- round(x$epi$prev.rgc[at], 3)
          prev.ugc <- round(x$epi$prev.ugc[at], 3)
          prev.rct <- round(x$epi$prev.rct[at], 3)
          prev.uct <- round(x$epi$prev.uct[at], 3)

          cat("\014")
          cat("\nEpidemic Simulation")
          cat("\n==============================")
          cat("\nTimestep: ", at, "/", nsteps, sep = "")
          cat("\nDay: ", at * time.unit, "/", nsteps * time.unit, sep = "")
          cat("\nYear: ", round((at * time.unit) / 365, 1), "/",
              round((nsteps * time.unit) / 365, 1), sep = "")
          cat("\n------------------------------")
          cat("\nTotal Pop Size:", x$epi$num[at])
          cat("\n------------------------------")
          cat("\nHIV Curr Incidence:", x$epi$incid[at])
          cat("\nHIV Cuml Incidence:", sum(x$epi$incid, na.rm = TRUE))
          cat("\nHIV Prevalence: ", x$epi$i.num[at], " (", prev, ")", sep = "")
          cat("\n------------------------------")
          cat("\nrGC Prevalence: ", prev.rgc, sep = "")
          cat("\nuGC Prevalence: ", prev.ugc, sep = "")
          cat("\nrCT Prevalence: ", prev.rct, sep = "")
          cat("\nuCT Prevalence: ", prev.uct, sep = "")
          cat("\n==============================")

        }
      }
    }


  }
}
