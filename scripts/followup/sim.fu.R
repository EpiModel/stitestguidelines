
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
stihrmodel <- as.character(Sys.getenv("STIHRMODEL"))
anncov <- as.numeric(Sys.getenv("ANNCOV"))
hrcov <- as.numeric(Sys.getenv("HRCOV"))
annint <- as.numeric(Sys.getenv("ANNINT"))
hrint <- as.numeric(Sys.getenv("HRINT"))
cov <- as.numeric(Sys.getenv("COV"))
prstiint <- as.numeric(Sys.getenv("PSTIINT"))
rc <- as.numeric(Sys.getenv("RC"))
probtx <- as.numeric(Sys.getenv("PROBTX"))
asymptx <- as.numeric(Sys.getenv("ASYMPTX"))

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,
                   ai.scale = 1.12, # 1.11889726, # was 1.13
                   
                   rsyph.tprob = 0.065,
                   usyph.tprob = 0.045,
                   
                   hiv.rsyph.rr = 3.33739617, 
                   hiv.usyph.rr = 1.85287705,
                   syph.rhiv.rr = 7.34516422,
                   syph.uhiv.rr = 3.64002970,
                   
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
                   ugc.tprob = 0.30904720, 
                   rct.tprob = 0.1907554, 
                   uct.tprob = 0.16394697,
                   
                   hiv.rgc.rr = 2.30,
                   hiv.ugc.rr = 1.30,
                   hiv.rct.rr = 2.30,
                   hiv.uct.rr = 1.30,
                   
                   syph.prim.sympt.prob.tx = 0.35, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.prim.asympt.prob.tx = 0.00,
                   syph.seco.sympt.prob.tx = 0.60, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.seco.asympt.prob.tx = 0.00,
                   syph.earlat.prob.tx = 0.15, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 0.90,
                   syph.tert.asympt.prob.tx = 0.00,

                   prep.sti.screen.int = prstiint,
                   prep.sti.prob.tx = probtx,
                   
                   gc.asympt.prob.tx = asymptx,
                   ct.asympt.prob.tx = asymptx,
                   
                   hivdx.syph.sympt.tx.rr = 1.45,
                   
                   stitest.elig.model = stihrmodel,
                   stianntest.coverage = anncov,
                   stihighrisktest.coverage = hrcov,
                   ept.coverage = 0,
                   
                   rcomp.prob = rc,
                   rcomp.adh.groups = 2:3,
                    
                   prep.start = 2601,
                   stitest.start = 2601,
                   
                   stitest.active.int = annint,
                   sti.highrisktest.int = hrint) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 2601,
                       nsteps = 3120,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = FALSE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = FALSE)
