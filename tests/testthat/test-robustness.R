# Robustness checks and predictor windows
# Converted from tests/test_robustness.R. Runs from the module folder regardless of the working directory.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))

testthat::test_that("Robustness checks and predictor windows", {
  withr::local_dir(root)
  ok <- function(cond, msg) testthat::expect_true(isTRUE(cond), info = msg)
  check <- function(name, cond) testthat::expect_true(isTRUE(cond), info = name)
  source("R/engine.R")
  d <- read.csv("data/california_prop99.csv", check.names = FALSE)
  Y <- as.matrix(d[, -1]); rownames(Y) <- d$Year; yrs <- as.integer(d$Year)
  pre <- which(yrs <= 1988); post <- which(yrs >= 1989)

  cat("Robustness checks (Prop 99)\n")
  main <- scm_fit(Y, "California", pre, post, starts = 5)
  att_main <- scm_att(main)
  cat(sprintf("  main ATT = %.3f\n", att_main))
  ok(att_main < -15, "main ATT is a large negative effect")

  t0 <- system.time(loo <- scm_loo(Y, "California", pre, post))[["elapsed"]]
  cat(sprintf("  leave-one-out time %.1fs\n", t0))
  ok(nrow(loo) == 38, "leave-one-out has 38 rows")
  ok(all(loo$att < 0), "ATT stays negative when any single donor is dropped")
  cat(sprintf("  LOO ATT range [%.2f, %.2f]\n", min(loo$att), max(loo$att)))

  tim <- scm_timing(Y, "California", 1989)
  ok(nrow(tim) == 4 && all(tim$att < 0), "timing sensitivity: 4 rows, all negative")
  pl <- scm_prelength(Y, "California", 1989)
  ok(nrow(pl) == 4, "pre-period length: 4 rows")
  it <- scm_intime(Y, "California", 1989)
  ok(nrow(it) == 3, "in-time placebo: 3 rows")
  cat("  in-time placebo ATTs:", round(it$att, 2), "\n")
  cat("  timing ATTs:", round(tim$att, 2), "\n")
  cat("  pre-length ATTs:", round(pl$att, 2), "\n")
  cat("Predictor summaries and window-specific predictors\n")
  # Tiny long panel: 2 units, 4 periods, predictor z with known values.
  long <- data.frame(unit = rep(c("a", "b"), each = 4), year = rep(1:4, 2),
                     z = c(1, 2, 3, 10,  5, 5, 5, 5))
  ok(isTRUE(all.equal(scm_predictor_matrix(long, "unit", "year", "z", c("a", "b"), 1:4, "mean"),
                      matrix(c(4, 5), 1, 2, dimnames = list("z", c("a", "b"))))),
     "mean over all pre periods")
  ok(isTRUE(all.equal(scm_predictor_matrix(long, "unit", "year", "z", c("a", "b"), 1:3, "last"),
                      matrix(c(3, 5), 1, 2, dimnames = list("z", c("a", "b"))))),
     "last pre value uses the latest period in the window")
  ok(isTRUE(all.equal(scm_predictor_matrix(long, "unit", "year", "z", c("a", "b"), 1:4, "median"),
                      matrix(c(2.5, 5), 1, 2, dimnames = list("z", c("a", "b"))))),
     "median over pre periods")
  ok(is.null(scm_predictor_matrix(long, "unit", "year", character(0), c("a", "b"), 1:4)),
     "no predictors gives NULL (outcome-only)")

  # Window-specific predictors: a predictor that jumps after the fake treatment year
  # must not leak into the placebo estimate. X_fun(pre) must ignore post-window rows.
  Xf <- function(pre_times) scm_predictor_matrix(long, "unit", "year", "z", c("a", "b"), pre_times, "mean")
  ok(isTRUE(all.equal(Xf(1:2)["z", "a"], 1.5)), "window predictor uses only the window's pre periods")
  ok(isTRUE(all.equal(Xf(1:4)["z", "a"], 4)), "full-window predictor differs from the short window")

  # Seed handling: scm_fit must not change the caller's RNG stream.
  set.seed(7); expected <- runif(1)
  set.seed(7); invisible(scm_fit(Y, "California", pre, post, starts = 2, seed = 3)); got <- runif(1)
  ok(identical(expected, got), "scm_fit restores the caller's RNG state")

  # Robustness checks run with predictors (X_fun path) and return finite ATTs.
  prop99_long <- data.frame(unit = rep(colnames(Y), each = nrow(Y)),
                            year = rep(as.integer(rownames(Y)), ncol(Y)),
                            cigsale = as.vector(Y), lnincome = 1 + 0.01 * seq_len(nrow(Y) * ncol(Y)))
  units_all <- colnames(Y)
  Xp <- function(pt) scm_predictor_matrix(prop99_long, "unit", "year", "lnincome", units_all, pt, "mean")
  tmp <- scm_timing(Y, "California", 1989, X_fun = Xp, starts = 1)
  ok(nrow(tmp) == 4 && all(is.finite(tmp$att)), "timing check runs with predictors")
  ptmp <- scm_prelength(Y, "California", 1989, drop = c(0, 4), X_fun = Xp, starts = 1)
  ok(nrow(ptmp) == 2 && all(is.finite(ptmp$att)), "pre-period length check runs with predictors")
  itmp <- scm_intime(Y, "California", 1989, X_fun = Xp, starts = 1)
  ok(nrow(itmp) == 3 && all(is.finite(itmp$att)), "in-time placebo runs with predictors")

  cat("\nAll robustness tests passed.\n")

})
