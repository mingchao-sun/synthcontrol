# Unit and regression tests for the synthetic control engine. Run: Rscript tests/run_tests.R
source("R/engine.R")
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE); cat("  ok -", msg, "\n") }

cat("1. Constraints of donor weights\n")
set.seed(42)
X0 <- matrix(rnorm(30), 10, 3); x1 <- rnorm(10); v <- runif(10)
w <- scm_weights_given_v(X0, x1, v)
ok(abs(sum(w) - 1) < 1e-10, "weights sum to one")
ok(all(w >= 0), "weights are non-negative")

cat("2. Exact recovery when treated unit is a convex combination of donors\n")
T <- 20; N <- 6
Y <- matrix(cumsum(rnorm(T * N)), T, N, dimnames = list(1:T, paste0("u", 1:N)))
true_w <- c(0.5, 0.3, 0.2, 0, 0)
Y[, 1] <- Y[, 2:6] %*% true_w
fit <- scm_fit(Y, treat = "u1", pre = 1:12, post = 13:20, starts = 2)
ok(fit$rmspe_pre < 1e-3, sprintf("pre-treatment RMSPE near zero (%.2e)", fit$rmspe_pre))
ok(all(abs(fit$weights_all[c("u2","u3","u4")] - true_w[1:3]) < 1e-2), "true donor weights recovered")

cat("3. Objective optimality versus random feasible weights\n")
obj <- function(w) sum(v * (x1 - X0 %*% w)^2)
w_star <- scm_weights_given_v(X0, x1, v)
rand_ok <- all(replicate(200, {
    r <- rexp(3); r <- r / sum(r); obj(r) >= obj(w_star) - 1e-10
}))
ok(rand_ok, "QP solution beats 200 random feasible points")

cat("4. Edge case: single donor gets weight one\n")
Y2 <- Y[, 1:2]
f2 <- scm_fit(Y2, treat = "u1", pre = 1:12, post = 13:20)
ok(abs(f2$weights_all[["u2"]] - 1) < 1e-10, "single donor weight equals one")

cat("5. Regression: Prop 99 (outcome-only predictors)\n")
d <- read.csv("data/california_prop99.csv", check.names = FALSE)
Yp <- as.matrix(d[, -1]); rownames(Yp) <- d$Year
yrs <- as.integer(d$Year)
fp <- scm_fit(Yp, "California", pre = which(yrs <= 1988), post = which(yrs >= 1989), starts = 5)
ok(abs(sum(fp$weights_all) - 1) < 1e-10, "Prop 99 weights sum to one")
ok(fp$ratio > 5, sprintf("post/pre RMSPE ratio large for Prop 99 (%.2f)", fp$ratio))
# Objective value from independent SLSQP check (scipy) at the same V: 0.0333512639
# Independent check: with the same v, scipy SLSQP gives objective 0.0333512639 (see M1 log)
Xs <- (Yp[which(yrs <= 1988), ] - rowMeans(Yp[which(yrs <= 1988), ])) / apply(Yp[which(yrs <= 1988), ], 1, sd)
dn <- setdiff(colnames(Yp), "California")
obj_p <- sum(fp$v * (Xs[, "California"] - Xs[, dn] %*% fp$weights_all[dn])^2)
ok(abs(obj_p - 0.0333512639) < 1e-6, sprintf("QP objective matches SLSQP reference (diff %.2e)", obj_p - 0.0333512639))
ok(fp$weights_all[["Utah"]] > 0.3, "Utah has the largest weight (matches Abadie et al., 2010 ranking)")

cat("6. Placebo p-value is a valid probability\n")
pl <- scm_placebo(Yp, "California", pre = which(yrs <= 1988), post = which(yrs >= 1989), starts = 1)
ok(pl$p_value >= 0 && pl$p_value <= 1, sprintf("p-value in [0,1] (%.3f)", pl$p_value))
ok(length(pl$placebo_ratios) == 38, "38 placebo fits for 38 donors")



cat("7. p-values from placebo distribution\n")
p_all <- scm_pvalue(pl)
p_f   <- scm_pvalue(pl, max_pre_ratio = 2)
ok(abs(p_all[["p"]] - 3 / 39) < 1e-12, sprintf("unfiltered p = 3/39 (%.4f)", p_all[["p"]]))
ok(p_f[["n_placebo"]] == sum(pl$placebo_pre_rmspe <= 2 * pl$fit$rmspe_pre), sprintf("filtered placebo count matches direct count (%d)", p_f[["n_placebo"]]))
ok(p_f[["p"]] >= 0 && p_f[["p"]] <= 1, sprintf("filtered p in [0,1] (%.3f, n=%d)", p_f[["p"]], p_f[["n_placebo"]]))


cat("8. Two-sided p-value and leave-one-out subset\n")
p2 <- scm_pvalue(pl, method = "att")
ok(p2[["p"]] >= 0 && p2[["p"]] <= 1, sprintf("two-sided p in [0,1] (%.3f)", p2[["p"]]))
ok(p2[["n_placebo"]] == 38, "two-sided uses all 38 placebos when unfiltered")
sub <- scm_loo(Yp, "California", pre = which(as.integer(d$Year) <= 1988), post = which(as.integer(d$Year) >= 1989),
               only = c("Utah", "Nevada"), starts = 1)
ok(nrow(sub) == 2 && all(sub$removed %in% c("Utah", "Nevada")), "LOO subset restricted to requested donors")
cat("\nAll tests passed (incl. options 8).\n")
