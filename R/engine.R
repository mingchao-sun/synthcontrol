# engine.R -----------------------------------------------------------------
# Pure-R synthetic control engine (Abadie, Diamond & Hainmueller 2010).
# Independent of jamovi / jmvcore so it can be unit-tested on its own.
#
# Conventions
#   Y      : numeric matrix, T x N (rows = time, cols = units). Column names = unit ids.
#   treat  : name or index of the treated unit (column of Y).
#   pre,post : integer row indices of pre- and post-treatment periods.
#   X      : optional K x N predictor matrix (rows = predictors, cols = units).
#            If NULL, the pre-treatment outcomes Y[pre, ] are used as predictors
#            ("outcome-only" specification).

require_quadprog <- function() {
    if (!requireNamespace("quadprog", quietly = TRUE))
        stop("Package 'quadprog' is required.", call. = FALSE)
}

# Donor weights for given predictor weights v (diagonal of V).
# min_w (x1 - X0 w)' diag(v) (x1 - X0 w)  s.t. w >= 0, sum(w) = 1
scm_weights_given_v <- function(X0, x1, v, ridge = 1e-10) {
    require_quadprog()
    J <- ncol(X0)
    Xv <- X0 * v                              # row-scaled by v
    D <- 2 * crossprod(X0, Xv)                # 2 X0' V X0
    d <- 2 * crossprod(Xv, x1)                # 2 X0' V x1
    D <- D + diag(ridge, J)                   # keep D positive definite
    Amat <- cbind(rep(1, J), diag(J))         # sum-to-one (equality) + nonnegativity
    bvec <- c(1, rep(0, J))
    sol <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = 1)
    w <- pmax(sol$solution, 0)
    w / sum(w)
}

# Pre-treatment MSPE criterion as a function of v (inner problem solved each time).
scm_loss_v <- function(lv, X0, x1, Y0pre, y1pre) {
    v <- lv
    w <- scm_weights_given_v(X0, x1, v)
    r <- y1pre - as.vector(Y0pre %*% w)
    sum(r^2)
}

# Full SCM fit for one treated unit.
scm_fit <- function(Y, treat, pre, post, X = NULL, starts = 5, seed = 1, maxit = 200) {
    require_quadprog()
    if (is.character(treat)) treat <- match(treat, colnames(Y))
    if (is.na(treat)) stop("Treated unit not found in data.", call. = FALSE)
    donors <- setdiff(seq_len(ncol(Y)), treat)

    if (is.null(X)) X <- Y[pre, , drop = FALSE]                # K x N, rows = pre periods
    X <- as.matrix(X)
    # standardise predictors across units so that a common v is meaningful
    s <- apply(X, 1, stats::sd); s[s == 0 | is.na(s)] <- 1
    Xs <- (X - rowMeans(X)) / s

    X0 <- Xs[, donors, drop = FALSE]
    x1 <- Xs[, treat]
    K <- nrow(Xs)
    Y0pre <- Y[pre, donors, drop = FALSE]
    y1pre <- Y[pre, treat]

    # outer optimisation over v in [0, 1]^K (multi-start; keep best).
    # The caller's RNG state is saved and restored so that set.seed() has no side effect.
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                    get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
        if (is.null(old_seed)) {
            if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
                rm(".Random.seed", envir = globalenv())
        } else assign(".Random.seed", old_seed, envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
    best <- NULL
    inits <- c(list(rep(1 / K, K)), lapply(seq_len(max(starts - 1, 0)), function(i) stats::runif(K)))
    for (v0 in inits) {
        opt <- stats::optim(v0, scm_loss_v, X0 = X0, x1 = x1, Y0pre = Y0pre, y1pre = y1pre,
                            method = "L-BFGS-B", lower = rep(0, K), upper = rep(1, K),
                            control = list(maxit = maxit))
        if (is.null(best) || opt$value < best$value) best <- opt
    }
    v <- best$par
    if (sum(v) == 0) v <- rep(1 / K, K)
    w <- scm_weights_given_v(X0, x1, v)

    names(w) <- colnames(Y)[donors]
    synth <- as.vector(Y[, donors, drop = FALSE] %*% w)
    gap <- Y[, treat] - synth
    names(synth) <- names(gap) <- rownames(Y)

    rmspe_pre  <- sqrt(mean(gap[pre]^2))
    rmspe_post <- sqrt(mean(gap[post]^2))

    balance <- data.frame(
        predictor   = if (is.null(rownames(X))) paste0("X", seq_len(K)) else rownames(X),
        treated     = X[, treat],
        synthetic   = as.vector(X[, donors, drop = FALSE] %*% w),
        row.names = NULL
    )
    list(
        weights   = w[w > 1e-6],                   # non-zero donor weights
        weights_all = w,
        v         = setNames(v, balance$predictor),
        synth     = synth,
        treated_y = Y[, treat],
        gap       = gap,
        rmspe_pre = rmspe_pre,
        rmspe_post = rmspe_post,
        ratio     = rmspe_post / rmspe_pre,
        balance   = balance,
        loss      = best$value,
        pre = pre, post = post, treat = treat
    )
}

# In-space placebo: re-run the fit treating each donor as if treated.
# The true treated unit is excluded from the placebo donor pool.
scm_placebo <- function(Y, treat, pre, post, X = NULL, fit_true = NULL, ...) {
    if (is.character(treat)) treat <- match(treat, colnames(Y))
    units <- setdiff(seq_len(ncol(Y)), treat)
    res <- lapply(units, function(u) {
        cols <- c(u, setdiff(seq_len(ncol(Y)), c(u, treat)))   # placebo treated first
        Yu <- Y[, cols, drop = FALSE]
        Xu <- if (is.null(X)) NULL else as.matrix(X)[, cols, drop = FALSE]
        scm_fit(Yu, treat = 1, pre = pre, post = post, X = Xu, ...)
    })
    names(res) <- colnames(Y)[units]
    ratios <- vapply(res, `[[`, numeric(1), "ratio")
    pre_rmspe <- vapply(res, `[[`, numeric(1), "rmspe_pre")
    # the treated unit is ranked against placebo ratios (including itself)
    if (is.null(fit_true)) fit_true <- scm_fit(Y, treat, pre, post, X = X, ...)
    all_ratios <- c(fit_true$ratio, ratios)
    p_value <- mean(all_ratios >= fit_true$ratio)
    list(
        fit = fit_true,
        placebo = res,
        placebo_ratios = ratios,
        placebo_pre_rmspe = pre_rmspe,
        p_value = p_value,
        gaps = do.call(cbind, lapply(res, `[[`, "gap"))
    )
}

# Permutation p-value from an scm_placebo() result.
# max_pre_ratio: keep only placebos whose pre-treatment RMSPE is at most
# max_pre_ratio times the treated unit's (Inf = no filter). The treated unit
# itself is always included in the ranking.
scm_pvalue <- function(pl, max_pre_ratio = Inf, method = c("ratio", "att")) {
    method <- match.arg(method)
    keep <- pl$placebo_pre_rmspe <= max_pre_ratio * pl$fit$rmspe_pre
    if (method == "ratio") {
        tr <- pl$fit$ratio
        r  <- pl$placebo_ratios[keep]
    } else {                       # two-sided: rank absolute average post-treatment gap
        tr <- abs(scm_att(pl$fit))
        r  <- abs(vapply(pl$placebo, scm_att, numeric(1)))[keep]
    }
    c(p = mean(c(tr, r) >= tr), n_placebo = length(r))
}

# Predictor matrix (K x N): each predictor summarised over the units' rows whose
# time is in pre_times. stat is "mean", "median" or "last". Returns NULL when
# there are no predictors (outcome-only specification).
scm_predictor_matrix <- function(df, unit_name, time_name, predictors, units,
                                 pre_times, stat = "mean") {
    if (length(predictors) == 0) return(NULL)
    summarise <- switch(stat,
        mean   = function(x, t) mean(x, na.rm = TRUE),
        median = function(x, t) stats::median(x, na.rm = TRUE),
        last   = function(x, t) x[which.max(t)],
        stop("Unknown predictor summary: ", stat, call. = FALSE))
    rows <- df[[time_name]] %in% pre_times
    pu <- as.character(df[[unit_name]][rows])
    pt <- df[[time_name]][rows]
    X <- vapply(predictors, function(p) {
        px <- df[[p]][rows]
        vapply(units, function(u) summarise(px[pu == u], pt[pu == u]), numeric(1))
    }, numeric(length(units)))
    X <- t(X)                                   # K x N
    dimnames(X) <- list(predictors, units)
    X
}

# ---- Robustness checks ----------------------------------------------------
# Average post-treatment gap (ATT) of a fit.
scm_att <- function(fit) mean(fit$gap[fit$post])

# Fit on a time window given as row indices of Y (pre first, then post).
# Rows outside the window are dropped, so only the window is used for weights.
# X_fun(pre_times) returns the predictor matrix computed from the window's
# pre-treatment periods only (NULL for outcome-only), so predictors never use
# data from the window's post-period.
scm_fit_window <- function(Y, treat, pre_idx, post_idx, X_fun = NULL, starts = 3, ...) {
    sub <- Y[c(pre_idx, post_idx), , drop = FALSE]
    yrs <- as.numeric(rownames(Y))
    X <- if (is.null(X_fun)) NULL else X_fun(yrs[pre_idx])
    scm_fit(sub, treat, pre = seq_along(pre_idx),
            post = length(pre_idx) + seq_along(post_idx), X = X, starts = starts, ...)
}

# Leave-one-out: drop each donor in turn and re-estimate.
scm_loo <- function(Y, treat, pre, post, X = NULL, only = NULL, starts = 3, ...) {
    if (is.character(treat)) treat <- match(treat, colnames(Y))
    tname <- colnames(Y)[treat]
    donors <- if (is.null(only)) setdiff(seq_len(ncol(Y)), treat) else match(only, colnames(Y))
    rows <- lapply(donors, function(j) {
        keep <- setdiff(seq_len(ncol(Y)), j)
        Xk <- if (is.null(X)) NULL else as.matrix(X)[, keep, drop = FALSE]
        f <- scm_fit(Y[, keep, drop = FALSE], tname, pre, post, X = Xk, starts = starts, ...)
        data.frame(removed = colnames(Y)[j], att = scm_att(f), ratio = f$ratio,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Treatment-timing sensitivity: move the treatment year by offsets.
scm_timing <- function(Y, treat, t0, offsets = c(-2, -1, 1, 2), X_fun = NULL, starts = 3, ...) {
    yrs <- as.numeric(rownames(Y))
    rows <- lapply(offsets, function(o) {
        t <- t0 + o
        f <- scm_fit_window(Y, treat, which(yrs < t), which(yrs >= t), X_fun = X_fun, starts = starts, ...)
        data.frame(check = "Treatment timing", scenario = sprintf("Treatment year %d", t),
                   x = t, att = scm_att(f), ratio = f$ratio, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Pre-period length sensitivity: drop the last d pre-treatment periods.
scm_prelength <- function(Y, treat, t0, drop = c(0, 2, 4, 6), X_fun = NULL, starts = 3, ...) {
    yrs <- as.numeric(rownames(Y))
    pre_all <- which(yrs < t0); post <- which(yrs >= t0)
    rows <- lapply(drop, function(d) {
        pre <- utils::head(pre_all, length(pre_all) - d)
        if (length(pre) < 2) return(NULL)
        f <- scm_fit_window(Y, treat, pre, post, X_fun = X_fun, starts = starts, ...)
        data.frame(check = "Pre-period length",
                   scenario = sprintf("Pre-period ends %d (%d years)", max(yrs[pre]), length(pre)),
                   x = max(yrs[pre]), att = scm_att(f), ratio = f$ratio, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# In-time placebo: pretend treatment happened earlier, using only pre-treatment data.
scm_intime <- function(Y, treat, t0, offsets = c(-6, -4, -2), X_fun = NULL, starts = 3, ...) {
    yrs <- as.numeric(rownames(Y))
    rows <- lapply(offsets, function(o) {
        f0 <- t0 + o
        f <- scm_fit_window(Y, treat, which(yrs < f0), which(yrs >= f0 & yrs < t0),
                            X_fun = X_fun, starts = starts, ...)
        data.frame(check = "In-time placebo", scenario = sprintf("Fake treatment year %d", f0),
                   x = f0, att = scm_att(f), ratio = f$ratio, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Parse a user-entered comma-separated list of numbers, e.g. "-2,-1,1,2".
scm_parse_list <- function(s, name = "list") {
    parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    v <- suppressWarnings(as.numeric(parts))
    if (length(v) == 0 || anyNA(v))
        stop(sprintf("'%s' must be a comma-separated list of numbers, e.g. -2,-1,1,2", name), call. = FALSE)
    v
}

# ---- panel construction and validation (used by the jamovi analysis and tests) --
# Converts a long data frame to the T x N outcome matrix and checks the design.
# Errors are signalled with stop(); the analysis converts them into jamovi messages.
scm_build_panel <- function(df, Y_name, unit_name, time_name, treated, treatmentTime,
                            exclude = character(0)) {
    # jamovi may read the time variable as nominal (factor/text).
    df[[time_name]] <- suppressWarnings(as.numeric(as.character(df[[time_name]])))
    if (anyNA(df[[time_name]]))
        stop("The time variable must be numeric (e.g. years such as 1970, 1971, ...).")

    if (length(exclude) > 0) {
        if (treated %in% exclude)
            stop("The treated unit cannot appear in the donor exclusion list.")
        bad <- setdiff(exclude, unique(as.character(df[[unit_name]])))
        if (length(bad) > 0)
            stop(paste0("Unknown unit(s) in the exclusion list: ", paste(bad, collapse = ", ")))
        df <- df[!(as.character(df[[unit_name]]) %in% exclude), , drop = FALSE]
    }

    units <- unique(as.character(df[[unit_name]]))
    if (length(units) < 3)
        stop("Need at least two donor units after exclusions.")
    if (!(treated %in% units))
        stop("The treated unit was not found in the unit variable.")
    times <- sort(unique(df[[time_name]]))

    # Duplicate unit-period rows would be silently overwritten below.
    if (anyDuplicated(paste(df[[unit_name]], df[[time_name]])))
        stop("Each unit must have exactly one row per time period (duplicate unit-time rows found).")

    Y <- matrix(NA_real_, nrow = length(times), ncol = length(units),
                dimnames = list(as.character(times), units))
    Y[cbind(match(df[[time_name]], times),
            match(as.character(df[[unit_name]]), units))] <- df[[Y_name]]

    if (anyNA(Y))
        stop("The panel must be balanced: every unit needs an outcome in every period.")

    pre  <- which(times <  treatmentTime)
    post <- which(times >= treatmentTime)
    if (length(pre) < 2 || length(post) < 1)
        stop("Need at least two pre-treatment and one post-treatment period.")

    list(Y = Y, pre = pre, post = post, units = units, times = times, df = df)
}
