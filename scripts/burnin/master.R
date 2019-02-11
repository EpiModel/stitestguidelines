
## build master.sh script ##

library("EpiModelHPC")


# Calibration/Testing -----------------------------------------------------

# vars <- list(gcadd = seq(0, 0.1, 0.005))
sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.sh",
              simno.start = 1000,
              ckpt = TRUE,
              nsims = 500,
              ncores = 28,
              walltime = "2:00:00",
              mem = "100G")

sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.ikt.sh",
              simno.start = 3000,
              ckpt = FALSE,
              nsims = 500,
              ncores = 16,
              walltime = "2:00:00",
              mem = "55G")
sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.ikt.sh",
              simno.start = 3001,
              append = TRUE,
              ckpt = TRUE,
              nsims = 500,
              ncores = 16,
              walltime = "2:00:00",
              mem = "55G")
sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.ikt.sh",
              simno.start = 3002,
              append = TRUE,
              ckpt = TRUE,
              nsims = 1000,
              ncores = 16,
              walltime = "2:00:00",
              mem = "55G")
sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.ikt.sh",
              simno.start = 3003,
              append = TRUE,
              ckpt = TRUE,
              nsims = 2000,
              ncores = 16,
              walltime = "2:00:00",
              mem = "55G")
sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.ikt.sh",
              simno.start = 3004,
              append = TRUE,
              ckpt = TRUE,
              nsims = 2000,
              ncores = 16,
              walltime = "2:00:00",
              mem = "55G")
sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.ikt.sh",
              simno.start = 3005,
              append = TRUE,
              ckpt = TRUE,
              nsims = 4000,
              ncores = 16,
              walltime = "2:00:00",
              mem = "55G")
