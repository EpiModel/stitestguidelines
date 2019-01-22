
## build master.sh script ##

library("EpiModelHPC")


# Calibration/Testing -----------------------------------------------------

sbatch_master(vars = NULL,
              master.file = "scripts/burnin/master.sh",
              build.runsim = FALSE,
              # param.sheet = "params.csv",
              # param.tag = "Calibrate",
              append = FALSE,
              simno.start = 100,
              ckpt = TRUE,
              nsims = 100,
              ncores = 16,
              mem = "55G")

