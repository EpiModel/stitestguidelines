
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- Sys.getenv("SIMNO")
jobno <- Sys.getenv("PBS_ARRAYID")
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

## Parameters
load("est/nwstats.rda")
#load("est/abc.syph.parms.rda")
#for (i in seq_along(mean.p)) {
#    assign(names(mean.p)[i], unname(mean.p[i]))
#}

param <- param_msm(nwstats = st,
                   ai.scale = 1.11, # 1.11889726, # was 1.13
                   prep.start = 5000,
                   stitest.start = 5000,
                   
                   hiv.rgc.rr = 2.5, #2.780673,
                   hiv.ugc.rr = 1.5, #1.732363,
                   hiv.rct.rr = 2.5, #2.780673,
                   hiv.uct.rr = 1.5, #1.732363,
                   hiv.rsyph.rr = 2.5, # 2.00892218, # 2.00793856, # functional at 2.00
                   hiv.usyph.rr = 1.5,
                   syph.hiv.rr = 1.7, #2.17038393, # 2.16933746, # functional at 2.40
                   
                   rsyph.tprob = 0.038, # 0.01950727,
                   usyph.tprob = 0.030, # 0.01950727,
                   
                   rgc.tprob = 0.40, #0.3928965, # 0.38353111, # was 0.357698 # functional at 0.40
                   ugc.tprob = 0.32, #0.24297633, # 0.25444490, # was 0.248095 # functional at 0.35
                   rct.tprob = 0.21, #0.29367628, # 0.31968155, # was 0.321597 # functional at 0.21
                   uct.tprob = 0.16, #0.25309465,# 0.23424104, # was 0.212965 # functional at 0.15
                   
                   prep.coverage = 0,
                   ept.coverage = 0,
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(nwstats = st, 
                 prev.B = 0.10, 
                 prev.W = 0.10,
                 prev.ugc = 0.013,
                 prev.rgc = 0.013,
                 prev.uct = 0.013,
                 prev.rct = 0.013,
                 prev.syph.B = 0.01,
                 prev.syph.W = 0.01)

control <- control_msm(simno = fsimno,
                       nsteps = 3120,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
