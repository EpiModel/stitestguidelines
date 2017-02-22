
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
                   ai.scale = 1.12, # 1.11889726, # was 1.13
                   
                   rsyph.tprob = 0.07,
                   usyph.tprob = 0.05,
                   
                   hiv.rsyph.rr = 2.62510586, 
                   hiv.usyph.rr = 1.64580886,
                   syph.rhiv.rr = 2.70,
                   syph.uhiv.rr = 1.70,
                   
                   syph.earlat.rr = 0.5, #2/3, 0
                   incu.syph.int = 27,
                   prim.syph.int = 45,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 120 - 45 - 27,
                   latelat.syph.int = 9*365,
                   latelatelat.syph.int = 30*365,
                   tert.syph.int = 30*365,
                   immune.syph.int = 5*365,
                   syph.tert.prog.prob = 0.15,
                   
                   rgc.tprob = 0.41, 
                   ugc.tprob = 0.31, 
                   rct.tprob = 0.20, 
                   uct.tprob = 0.15,
                   
                   hiv.rgc.rr = 2.62510586, #2.780673,
                   hiv.ugc.rr = 1.64580886, #1.732363,
                   hiv.rct.rr = 2.62510586, #2.780673,
                   hiv.uct.rr = 1.64580886, #1.732363,
                   
                   # adjust prim and seco from 0.1385 each
                   stage.syph.B.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   stage.syph.W.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   
                   syph.prim.sympt.prob.tx = 0.3771634, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.prim.asympt.prob.tx = 0.00,
                   syph.seco.sympt.prob.tx = 0.69365575, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.seco.asympt.prob.tx = 0.00,
                   syph.earlat.prob.tx = 0.16156070, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.prob.tx = 0.12667613,
                   syph.tert.sympt.prob.tx = 0.90,
                   syph.tert.asympt.prob.tx = 0.00,
                   
                   prep.coverage = 0,
                   ept.coverage = 0,
                   
                   prep.start = 5000,
                   stitest.start = 5000,
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(nwstats = st, 
                 prev.B = 0.10, 
                 prev.W = 0.10,
                 prev.ugc = 0.015,
                 prev.rgc = 0.015,
                 prev.uct = 0.015,
                 prev.rct = 0.015,
                 prev.syph.B = 0.02,
                 prev.syph.W = 0.02)

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = FALSE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")
