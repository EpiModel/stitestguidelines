
## Process burn-in
library("EpiModelHPC")

# Examine output
system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")

list.files("data/")

load("data/sim.n100.rda")
plot(sim, y = "i.prev", ylim = c(0.2, 0.3), qnts = 0.5)
abline(h = 0.26)

df <- as.data.frame(sim)
round(mean(tail(df$i.prev, 100)), 3)

# Set targets for multivariate Euclidean distance calculation
# rect.prev, ureth.prev, gc.incid, ct.incid, hiv.incid, hiv.prev
targets <- c(0.135, 0.046, 23.2, 26.8, 3.8, 0.26)

# Merge sim files
sim <- merge_simfiles(simno = 2200, indir = "data/", ftype = "max")

# Create function for selecting sim closest to target
mean_sim <- function(sim, targets, return.sim = TRUE, ...) {
    
    if (class(sim) != "netsim") {
        stop("sim must be of class netsim", call. = FALSE)
    }
    nsims <- sim$control$nsims
    if (missing(nsims)) {
        stop("Specify sims as a vector of simulations or \"mean\" ",
             call. = FALSE)
    }
    
    # Initialize distance vector
    dist <- c()
    calib <- c()
    targets <- as.data.frame(targets)
    
    # Obtain statistics and perform multivariable Euclidean distance calculation in a for-loop
    for (i in 1:nsims){
        
        # Create data frame to draw statistics from
        df <- as.data.frame(x=sim, out="vals", sim = i)
        
        # Create a vector of statistics
        calib[i] <- as.data.frame(c(mean(tail(df$prev.rgcct, 250)), mean(tail(df$prev.ugcct, 250)), mean(tail(df$ir100.gc, 250)),
                                    mean(tail(df$ir100.ct, 250)), mean(tail(df$ir100, 250)), mean(tail(df$i.prev, 250))))
        
        # Iteratively calculate distance
        dist[i] <- sqrt(sum((calib[i] - targets)^2))
    }
    
    # Which sim minimizes distance
    meansim <- which.min(dist)
    delsim <- setdiff(1:nsims, meansim)
    out <- sim
    if (length(delsim) > 0) {
        for (i in seq_along(out$epi)) {
            out$epi[[i]] <- out$epi[[i]][, -delsim, drop = FALSE]
        }
        if (!is.null(out$network)) {
            out$network[delsim] <- NULL
        }
        if (!is.null(out$stats$nwstats)) {
            out$stats$nwstats[delsim] <- NULL
        }
        if (!is.null(out$stats$transmat)) {
            out$stats$transmat[delsim] <- NULL
        }
        if (!is.null(out$control$save.other)) {
            oname <- out$control$save.other
            for (i in seq_along(oname)) {
                out[[oname[i]]][delsim] <- NULL
            }
        }
    }
    out$control$nsims <- length(meansim)
    
    # Print table of targets and statistics from the selected sim
    diff <- calib[meansim]-targets
    cat(meansim)
    print(data.frame(calib[meansim], targets, diff, dist[meansim]))
    
    # Returning simulation
    if(return.sim == TRUE) {
        sim <<- get_sims(sim, meansim)
    }
}

# Run function
mean_sim(sim, targets)


# Save burn-in file for FU sims
sim <- merge_simfiles(simno = 100, indir = "data/", ftype = "max")
sim <- get_sims(sim, sims = "mean", var = "i.prev")
tail(as.data.frame(sim)$i.prev)

save(sim, file = "est/stimod.burnin.rda")
system("scp est/stimod.burnin.rda hyak:/gscratch/csde/sjenness/sti/est/")
