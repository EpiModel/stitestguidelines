
## build master.sh script ##

# base model
vars <- list(COV = 0,
             ADR = 0,
             RC1 = 0,
             RC2 = FALSE,
             RC3 = FALSE,
             RC4 = FALSE)
qsub_master(simno.start = 2000,
            nsubjobs = 8,
            backfill = FALSE,
            vars = vars,
            append = FALSE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

# RC * high-adherent only interaction model
vars <- list(COV = 0.4,
             ADR = 0.019,
             RC1 = seq(0, 1, 0.05),
             RC2 = c(FALSE, TRUE),
             RC3 = FALSE,
             RC4 = FALSE)
qsub_master(simno.start = 2100,
            nsubjobs = 8,
            backfill = 15,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

# RC when with main partners only
vars <- list(COV = 0.4,
             ADR = 0.019,
             RC1 = seq(0, 1, 0.05),
             RC2 = FALSE,
             RC3 = TRUE,
             RC4 = FALSE)
qsub_master(simno.start = 2200,
            nsubjobs = 8,
            backfill = 15,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

# RC when with known serodiscordant only
vars <- list(COV = 0.4,
             ADR = 0.019,
             RC1 = seq(0, 1, 0.05),
             RC2 = FALSE,
             RC3 = FALSE,
             RC4 = TRUE)
qsub_master(simno.start = 2300,
            nsubjobs = 8,
            backfill = 15,
            vars = vars,
            append = TRUE,
            runsimfile = "runsim.fu.sh",
            outfile = "scripts/followup/master.fu.sh")

