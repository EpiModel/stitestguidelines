load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   # STI acquisition (ABC)
                    rgc.tprob = 0.57, #ABC Range: 0.54 - 0.60
                    ugc.tprob = 0.50, #ABC Range: 0.47 - 0.53
                    rct.tprob = 0.245,#ABC Range: 0.215 - 0.275
                    uct.tprob = 0.205,#ABC Range: 0.175 - 0.235
                    syph.tprob = 0.26,

                    # Possible ABC (if not computationally burdensome)
                    rgc.asympt.int = 21*7, # ABC Range: 15 -25 weeks (e.g. 15*7 - 25.7)
                    ugc.asympt.int = 21*7, # ABC Range: 15 -25 weeks (e.g. 15*7 - 25.7)

                    # Reviewer requested edits
                    rgc.sympt.prob = 0.16, # Beck
                    ugc.sympt.prob = 0.80, # Beck value - 0.10 (reviewer)
                    rct.sympt.prob = 0.14, # Beck
                    uct.sympt.prob = 0.48, # Beck value - 0.10 (reviewer)

                    #### Other parameters

                    # Correlation
                    sti.correlation.time = 12,

                    # HIV acquisition
                    hiv.rgc.rr = 1.97, #1.75,
                    hiv.ugc.rr = 1.48, #1.27,
                    hiv.rct.rr = 1.97, #1.75,
                    hiv.uct.rr = 1.48, #1.27,
                    hiv.syph.rr = 1.64,

                    # Syphilis probabilities
                    syph.incub.sympt.prob = 0,
                    syph.prim.sympt.prob = 0.82,
                    syph.seco.sympt.prob = 0.90,
                    syph.earlat.sympt.prob = 0,
                    syph.latelat.sympt.prob = 0,
                    syph.tert.sympt.prob = 1.0,

                    syph.prim.sympt.prob.tx = 0.85,
                    syph.seco.sympt.prob.tx = 0.85,
                    syph.earlat.sympt.prob.tx = 0.10,
                    syph.latelat.sympt.prob.tx = 0.10,
                    syph.tert.sympt.prob.tx = 1.0,

                    # Intervention settings
                    ept.coverage = 0.0,
                    stianntest.gc.hivneg.coverage = 0.44,
                    stianntest.ct.hivneg.coverage = 0.44,
                    stianntest.syph.hivneg.coverage = 0.45,
                    stihighrisktest.gc.hivneg.coverage = 0,
                    stihighrisktest.ct.hivneg.coverage = 0,
                    stihighrisktest.syph.hivneg.coverage = 0,
                    stianntest.gc.hivpos.coverage = 0.61,
                    stianntest.ct.hivpos.coverage = 0.61,
                    stianntest.syph.hivpos.coverage = 0.67,
                    stihighrisktest.gc.hivpos.coverage = 0,
                    stihighrisktest.ct.hivpos.coverage = 0,
                    stihighrisktest.syph.hivpos.coverage = 0,

                    prep.start = 7000,
                    stitest.start = 1,
                    ept.start = 5201,

                    stitest.active.int = 364,
                    sti.highrisktest.int = 182,
                    ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.ugc = 0.0015,
                 prev.rgc = 0.0015,
                 prev.uct = 0.0015,
                 prev.rct = 0.0015, # 0.03
                 prev.syph.B = 0, # 0.03
                 prev.syph.W = 0) # 0.03

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 8, ncores = 8,
                       verbose = FALSE)

netsim_hpc("est/fit.rda", param, init, control,
           save.min = FALSE, save.max = TRUE)
