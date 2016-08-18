
## build master.sh script ##

library("EpiModelHPC")

vars <- list(COV = 0,
             PSTIINT = seq(7, 7*26, 14),
             RC = 0,
             ASYMPT = 0)
qsub_master(simno.start = 1000,
            nsubjobs = 2,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")


vars <- list(COV = c(0, 0.5, 0.9),
             PSTIINT = 182,
             RC = c(0, 0.5, 1))
qsub_master(simno.start = 1000,
            nsubjobs = 4,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")
