
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
# stihrmodel <- as.character(Sys.getenv("STIHRMODEL"))
# anncov <- as.numeric(Sys.getenv("ANNCOV"))
# hrcov <- as.numeric(Sys.getenv("HRCOV"))
# annint <- as.numeric(Sys.getenv("ANNINT"))
# hrint <- as.numeric(Sys.getenv("HRINT"))
# cov <- as.numeric(Sys.getenv("COV"))
# stiasymptx <- as.numeric(Sys.getenv("STIASYMPTX"))

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   ai.scale = 1.11,
                   
                   rsyph.tprob = 0.04668348,
                   usyph.tprob = 0.03598792,
                   
                   hiv.rsyph.rr = 2.98876572, 
                   hiv.usyph.rr = 1.7456618,
                   syph.rhiv.rr = 6.54189295,
                   syph.uhiv.rr = 5.09641658,
                   
                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.15 / (52 * 7 * 20),
                   
                   rgc.tprob = 0.4133300,
                   ugc.tprob = 0.31404720,
                   rct.tprob = 0.1907554,
                   uct.tprob = 0.16394697,
                   
                   hiv.rgc.rr = 2.35,
                   hiv.ugc.rr = 1.35,
                   hiv.rct.rr = 2.35,
                   hiv.uct.rr = 1.35,
                   
                   syph.prim.sympt.prob.tx = 0.35, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.seco.sympt.prob.tx = 0.60, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.earlat.sympt.prob.tx = 0.15, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 0.90,
                   
                   # syph.prim.asympt.prob.tx = stiasymptx,
                   # syph.seco.asympt.prob.tx = stiasymptx,
                   # syph.earlat.asympt.prob.tx = stiasymptx, 
                   # syph.latelat.asympt.prob.tx = stiasymptx,
                   # syph.tert.asympt.prob.tx = stiasymptx,
                   # gc.asympt.prob.tx = stiasymptx,
                   # ct.asympt.prob.tx = stiasymptx,
                  
                   hivdx.syph.sympt.tx.rr = 1.45,
                   
                   # stitest.elig.model = stihrmodel,
                   # stianntest.coverage = anncov,
                   # stihighrisktest.coverage = hrcov,
                   # prep.coverage = cov,
                   ept.coverage = 0,
                    
                   prep.start = 5500,
                   stitest.start = 5500)#,
                   
                   # stitest.active.int = annint,
                   # sti.highrisktest.int = hrint) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 5200,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE)
