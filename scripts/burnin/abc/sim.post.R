
suppressMessages(library("EpiABC"))

out <- get_posterior(wave = 27, input = "scripts/burnin/abc/data/")
summary(out, digits = 3)

boxplot(out, type = "stats", stats = 1:15)
boxplot(out, type = "stats", stats = 16:30)
boxplot(out, type = "stats", stats = 31:45)

plot(out, type = "stats", stats = 1:15, stats.positive = FALSE)
plot(out, type = "stats", stats = 16:30, stats.positive = FALSE)
plot(out, type = "stats", stats = 31:45, stats.positive = FALSE)

plot(out, type = "param")
boxplot(out, type = "param")

