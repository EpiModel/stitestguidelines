
## build master.sh script ##

library(EpiModelHPC)

vars <- list(COV = seq(0, 1, 0.1))
qsub_master(simno.start = 1000,
            nsubjobs = 7,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")


vars <- list(COV = c(0, 0.5, 0.9),
             PSTIINT = c(91, 182))
qsub_master(simno.start = 1000,
            nsubjobs = 4,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")
