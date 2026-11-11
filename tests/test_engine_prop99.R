# Diagnostic script; run from the package root: Rscript tests/test_engine_prop99.R
source("R/engine.R")

d <- read.csv("data/california_prop99.csv", check.names = FALSE)
Y <- as.matrix(d[, -1]); rownames(Y) <- d$Year
years <- as.integer(d$Year)
pre  <- which(years <= 1988)
post <- which(years >= 1989)

cat("T =", nrow(Y), " N =", ncol(Y), " pre =", length(pre), " post =", length(post), "\n")

t0 <- proc.time()
fit <- scm_fit(Y, treat = "California", pre = pre, post = post, starts = 5)
cat("fit time (s):", round((proc.time() - t0)[["elapsed"]], 2), "\n")

cat("\n-- weights (non-zero) --\n"); print(round(sort(fit$weights, decreasing = TRUE), 4))
cat("sum w =", sum(fit$weights_all), " min w =", min(fit$weights_all), "\n")
cat("RMSPE pre =", round(fit$rmspe_pre, 3), " post =", round(fit$rmspe_post, 3),
    " ratio =", round(fit$ratio, 2), "\n")


cat("\n-- placebo (in-space) --\n")
t1 <- proc.time()
pl <- scm_placebo(Y, treat = "California", pre = pre, post = post, starts = 2)
cat("placebo time (s):", round((proc.time() - t1)[["elapsed"]], 2), "\n")
cat("p-value (RMSPE ratio rank) =", round(pl$p_value, 4), "\n")
print(round(sort(pl$placebo_ratios, decreasing = TRUE)[1:5], 2))
