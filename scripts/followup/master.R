
## build master.sh script ##

library("EpiModelHPC")

system("scp scripts/followup/*.[Rs]* hyak:/gscratch/csde/sjenness/sti/")

vars <- list(COV = 0.4,
             PSTIINT = c(7, 14, 28, 56, 91, 182, 364),
             RC = 0)
qsub_master(simno.start = 1000,
            nsubjobs = 7,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")
