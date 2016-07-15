
library("EpiModelHIV")


load("scripts/burnin/simDataAll.rda")

# rgc, ugc, rct, uct
targets <- c(0.102, 0.111, 0.141, 0.084)
targets <- c(0.083, 0.015, 0.118, 0.027)

diff.rgc <- abs(sim$rgc.prev - targets[1])
diff.ugc <- abs(sim$ugc.prev - targets[2])
diff.rct <- abs(sim$rct.prev - targets[3])
diff.uct <- abs(sim$uct.prev - targets[4])

th <- 0.03

choice <- which(diff.rgc <= th & diff.ugc <= th & diff.rct <= th & diff.uct <= th)
simChosen <- sim[choice, ]
simChosen

cbind(colMeans(simChosen))

hist(sim$rgc.prev)
hist(sim$ugc.prev)
hist(sim$rct.prev)
hist(sim$uct.prev)

# For V1 --------------------------------------------------------------


rejection <- function(sim,
                      target.stat = c(0.26, 0),
                      threshold = 0.001) {

  diff1 <- abs(sim$stat.mean - target.stat[1])

  accepted <- which(diff1 <= threshold)

  post <- sim[accepted, ]
  return(post)
}




post <- rejection(fsim)
post

hist(post$ai.scale)
plot(density(post$ai.scale))
plot(density(post$stat.mean))


# For V2 --------------------------------------------------------------

# this plot show how many simulations under the threshold we need before stabilization
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(density(simChosen$ai.scale[1:10]), ylim = c(0, 35),
     col = "blue", xlim = c(1.2, 1.35))
for (i in 11:126) {
  lines(density(simChosen$ai.scale[1:i]),
        col = ifelse(i < 126, grey(1 - (i/nrow(simChosen))), "red"))
  Sys.sleep(0.25)
}
abline(v = mean(simChosen$ai.scale))
abline(v = median(simChosen$ai.scale), col = "red")

a <- sapply(1:nrow(simChosen), function(x) mean(simChosen$ai.scale[1:x]))
diff(a)
plot(diff(a), type = "l")
abline(h = 0, lty = 2)


# this plot shows what level of threshold we might choose if using a fixed distance
vec.thresh <- rev(seq(0.0001, 0.01, 0.0001))
for (i in 1:length(vec.thresh)) {
  simThresh <- rejection(sim, threshold = vec.thresh[i])
  if (i == 1) {
    plot(density(simThresh$ai.scale), ylim = c(0, 35),
         col = "blue", xlim = c(1.2, 1.35))
    abline(v = mean(simThresh$ai.scale), col = "blue")
  } else {
    col.choice <- ifelse(i < length(vec.thresh), grey(1 - (i/length(vec.thresh))), "red")
    lines(density(simThresh$ai.scale),
          col = col.choice)
    abline(v = mean(simThresh$ai.scale), col = col.choice)
    Sys.sleep(0.1)
  }
}


# this plot shows what level of threshold we might choose if using a quantile
rejection2 <- function(sim, target.stat = 0.26, threshold = 0.01) {
  edist <- sapply(1:nrow(sim), function(x) sqrt(sum((target.stat - sim$stat.mean[x])^2)))
  edist.quant <- quantile(edist, threshold)
  accepted <- which(edist <= edist.quant)
  post <- sim[accepted, ]
  return(post)
}

vec.thresh <- rev(seq(0.001, 0.1, 0.001))
for (i in 1:length(vec.thresh)) {
  simThresh <- rejection2(sim, threshold = vec.thresh[i])
  if (i == 1) {
    plot(density(simThresh$ai.scale), ylim = c(0, 35),
         col = "blue", xlim = c(1.2, 1.35))
    abline(v = mean(simThresh$ai.scale), col = "blue")
  } else {
    col.choice <- ifelse(i < length(vec.thresh), grey(1 - (i/length(vec.thresh))), "red")
    lines(density(simThresh$ai.scale),
          col = col.choice)
    abline(v = mean(simThresh$ai.scale), col = col.choice)
    Sys.sleep(0.1)
  }
}

mean(rejection2(sim, threshold = 0.1)$ai.scale)
mean(rejection2(sim, threshold = 0.05)$ai.scale)
mean(rejection2(sim, threshold = 0.01)$ai.scale)
mean(rejection2(sim, threshold = 0.005)$ai.scale)
mean(rejection2(sim, threshold = 0.001)$ai.scale)

a <- sapply(vec.thresh, function(x) mean(rejection2(sim, threshold = x)$ai.scale))
da <- diff(a)
plot(a, type = 'l')
plot(a, type = "l", ylim = c(1.2, 1.3))
plot(da, type = "l")
