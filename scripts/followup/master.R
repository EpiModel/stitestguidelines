
## build master.sh script ##

# interaction model
vars <- list(ADR = seq(-0.5, 0.3, 0.1),
             RC1 = seq(0, 1, 0.1),
             RC2 = c(FALSE, TRUE),
             RC3 = FALSE,
             RC4 = FALSE)
qsub_master(simno.start = 2000,
            nsubjobs = 8,
            backfill = TRUE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.sh")
