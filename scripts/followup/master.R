
## build master.sh script ##

library("EpiModelHPC")

system("scp scripts/followup/*.[Rs]* hyak:/gscratch/csde/sjenness/sti/")

vars <- list(COV = 0,
             PSTIINT = 182,
             RC = 0)
qsub_master(simno.start = 1000,
            nsubjobs = 16,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

vars <- list(COV = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
             PSTIINT = c(91, 182, 364),
             RC = c(0, 0.25, 0.5, 0.75, 1))
qsub_master(simno.start = "auto",
            nsubjobs = 16,
            backfill = TRUE,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")
