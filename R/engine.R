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

# Number formatting and row labels shared by the engine and the jamovi tables,
# so that table rows built in .init() match the rows filled in .run().
scm_fmt <- function(x) format(x, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)
scm_lab_timing   <- function(t) sprintf("Treatment year %s", scm_fmt(t))
scm_lab_prelength <- function(d) sprintf("Drop last %s pre-treatment period(s)", scm_fmt(d))
scm_lab_intime   <- function(t) sprintf("Fake treatment year %s", scm_fmt(t))

# In-space placebo: re-run the fit treating each donor as if treated.
# The true treated unit is excluded from the placebo donor pool.
# checkpoint() is called before each placebo fit so jamovi can stay responsive.
scm_placebo <- function(Y, treat, pre, post, X = NULL, fit_true = NULL,
                        checkpoint = function() NULL, ...) {
    if (is.character(treat)) treat <- match(treat, colnames(Y))
    if (is.null(fit_true)) fit_true <- scm_fit(Y, treat, pre, post, X = X, ...)
    units <- setdiff(seq_len(ncol(Y)), treat)
    res <- lapply(units, function(u) {
        checkpoint()
        cols <- c(u, setdiff(seq_len(ncol(Y)), c(u, treat)))   # placebo unit first
        Yu <- Y[, cols, drop = FALSE]
        Xu <- if (is.null(X)) NULL else as.matrix(X)[, cols, drop = FALSE]
        scm_fit(Yu, treat = 1, pre = pre, post = post, X = Xu, ...)
    })
    names(res) <- colnames(Y)[units]
    ratios <- vapply(res, `[[`, numeric(1), "ratio")
    list(
        fit = fit_true,
        p_value = mean(c(fit_true$ratio, ratios) >= fit_true$ratio),   # unfiltered, one-sided
        placebo = res,
        placebo_ratios = vapply(res, `[[`, numeric(1), "ratio"),
        placebo_pre_rmspe = vapply(res, `[[`, numeric(1), "rmspe_pre"),
        placebo_att = vapply(res, scm_att, numeric(1)),
        gaps = do.call(cbind, lapply(res, `[[`, "gap"))
    )
}

# Only the numbers the inference and the gap figure need, so the jamovi state
# stays small (the full placebo fits are not stored).
scm_placebo_slim <- function(pl) {
    list(ratio = pl$fit$ratio, rmspe_pre = pl$fit$rmspe_pre, att = scm_att(pl$fit),
         placebo_ratios = pl$placebo_ratios, placebo_pre_rmspe = pl$placebo_pre_rmspe,
         placebo_att = pl$placebo_att, gaps = pl$gaps)
}

# Permutation p-value. Accepts scm_placebo() output or its slim version.
# max_pre_ratio: keep only placebos whose pre-treatment RMSPE is at most
# max_pre_ratio times the treated unit's (Inf = no filter).
# The treated unit itself is always included in the ranking.
scm_pvalue <- function(pl, max_pre_ratio = Inf, method = c("ratio", "att")) {
    method <- match.arg(method)
    if (is.null(pl$ratio)) pl <- scm_placebo_slim(pl)
    keep <- pl$placebo_pre_rmspe <= max_pre_ratio * pl$rmspe_pre
    if (method == "ratio") {
        tr <- pl$ratio
        r  <- pl$placebo_ratios[keep]
    } else {                        # two-sided: absolute average post-treatment gap
        tr <- abs(scm_att_value(pl))
        r  <- abs(pl$placebo_att)[keep]
    }
    c(p = mean(c(tr, r) >= tr), n_placebo = length(r))
}
scm_att_value <- function(pl) if (!is.null(pl$att)) pl$att else scm_att(pl$fit)

# Predictor matrix (K x N): each predictor summarised over the units' rows whose
# time is in pre_times. stat is "mean", "median" or "last". Returns NULL when
# there are no predictors. Stops with a message that names the predictor and
# the units when a value cannot be computed (for example a unit with no
# non-missing value in the pre-treatment periods).
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
        vapply(units, function(u) {
            sel <- pu == u
            if (!any(sel)) NA_real_ else summarise(px[sel], pt[sel])
        }, numeric(1))
    }, numeric(length(units)))
    X <- t(X)                                   # K x N
    dimnames(X) <- list(predictors, units)
    bad <- which(!is.finite(X), arr.ind = TRUE)
    if (nrow(bad) > 0) {
        msgs <- vapply(unique(bad[, 1]), function(k) {
            us <- units[unique(bad[bad[, 1] == k, 2])]
            shown <- if (length(us) > 5) paste(c(us[1:5], "..."), collapse = ", ")
                     else paste(us, collapse = ", ")
            sprintf("Predictor '%s' has no usable values in the pre-treatment periods for: %s.",
                    predictors[k], shown)
        }, character(1))
        stop(paste(msgs, collapse = " "), call. = FALSE)
    }
    X
}

# ---- Robustness checks ----------------------------------------------------
# Average post-treatment gap (ATT) of a fit.
scm_att <- function(fit) mean(fit$gap[fit$post])

# Fit on a time window given as row indices of Y. Rows outside the window are
# dropped, so only the window is used for the weights. X_fun(pre_times) returns
# the predictor matrix computed from the window's pre-treatment periods only.
scm_fit_window <- function(Y, treat, pre_idx, post_idx, X_fun = NULL, starts = 3, ...) {
    sub <- Y[c(pre_idx, post_idx), , drop = FALSE]
    yrs <- as.numeric(rownames(Y))
    X <- if (is.null(X_fun)) NULL else X_fun(yrs[pre_idx])
    scm_fit(sub, treat, pre = seq_along(pre_idx),
            post = length(pre_idx) + seq_along(post_idx), X = X, starts = starts, ...)
}

# Leave-one-out: drop each donor in turn and re-estimate.
scm_loo <- function(Y, treat, pre, post, X = NULL, only = NULL, starts = 3,
                    checkpoint = function() NULL, ...) {
    if (is.character(treat)) treat <- match(treat, colnames(Y))
    tname <- colnames(Y)[treat]
    donors <- if (is.null(only)) setdiff(seq_len(ncol(Y)), treat) else match(only, colnames(Y))
    rows <- lapply(donors, function(j) {
        checkpoint()
        keep <- setdiff(seq_len(ncol(Y)), j)
        Xk <- if (is.null(X)) NULL else as.matrix(X)[, keep, drop = FALSE]
        f <- scm_fit(Y[, keep, drop = FALSE], tname, pre, post, X = Xk, starts = starts, ...)
        data.frame(removed = colnames(Y)[j], att = scm_att(f), ratio = f$ratio,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Treatment-timing sensitivity: move the treatment year by offsets.
scm_timing <- function(Y, treat, t0, offsets = c(-2, -1, 1, 2), X_fun = NULL, starts = 3,
                       checkpoint = function() NULL, ...) {
    yrs <- as.numeric(rownames(Y))
    rows <- lapply(offsets, function(o) {
        checkpoint()
        t <- t0 + o
        f <- scm_fit_window(Y, treat, which(yrs < t), which(yrs >= t),
                            X_fun = X_fun, starts = starts, ...)
        data.frame(scenario = scm_lab_timing(t), x = t, att = scm_att(f), ratio = f$ratio,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Pre-period length sensitivity: drop the last d pre-treatment periods.
scm_prelength <- function(Y, treat, t0, drop = c(0, 2, 4, 6), X_fun = NULL, starts = 3,
                          checkpoint = function() NULL, ...) {
    yrs <- as.numeric(rownames(Y))
    pre_all <- which(yrs < t0); post <- which(yrs >= t0)
    rows <- lapply(drop, function(d) {
        checkpoint()
        pre <- utils::head(pre_all, length(pre_all) - d)
        f <- scm_fit_window(Y, treat, pre, post, X_fun = X_fun, starts = starts, ...)
        data.frame(scenario = scm_lab_prelength(d), x = max(yrs[pre]), att = scm_att(f),
                   ratio = f$ratio, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# In-time placebo: pretend treatment happened earlier, using only pre-treatment data.
scm_intime <- function(Y, treat, t0, offsets = c(-6, -4, -2), X_fun = NULL, starts = 3,
                       checkpoint = function() NULL, ...) {
    yrs <- as.numeric(rownames(Y))
    rows <- lapply(offsets, function(o) {
        checkpoint()
        f0 <- t0 + o
        f <- scm_fit_window(Y, treat, which(yrs < f0), which(yrs >= f0 & yrs < t0),
                            X_fun = X_fun, starts = starts, ...)
        data.frame(scenario = scm_lab_intime(f0), x = f0, att = scm_att(f), ratio = f$ratio,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Check every requested robustness scenario against the data before estimation.
# Errors name the offending value.
scm_check_scenarios <- function(times, t0, timing = NULL, drops = NULL, intime = NULL) {
    for (o in timing) {
        t <- t0 + o
        n_pre <- sum(times < t); n_post <- sum(times >= t)
        if (n_pre < 2 || n_post < 1)
            stop(sprintf(paste("Treatment timing offset %s moves the treatment to %s, which leaves",
                               "%d earlier and %d later periods. Each scenario needs at least two",
                               "earlier and one later period."),
                         scm_fmt(o), scm_fmt(t), n_pre, n_post), call. = FALSE)
    }
    if (length(drops) > 0) {
        if (any(drops < 0) || any(drops != round(drops)))
            stop("Pre-period drops must be non-negative whole numbers of periods.", call. = FALSE)
        n_pre_all <- sum(times < t0)
        if (any(n_pre_all - drops < 2))
            stop(sprintf(paste("Pre-period drop %s leaves fewer than two pre-treatment periods",
                               "(there are %d)."), scm_fmt(max(drops)), n_pre_all), call. = FALSE)
    }
    for (o in intime) {
        f0 <- t0 + o
        n_pre <- sum(times < f0); n_post <- sum(times >= f0 & times < t0)
        if (n_post < 1)
            stop(sprintf(paste("In-time placebo offset %s leaves no periods between the placebo date",
                               "and treatment; use negative offsets such as -6,-4,-2."), scm_fmt(o)),
                 call. = FALSE)
        if (n_pre < 2)
            stop(sprintf("In-time placebo offset %s leaves fewer than two pre-treatment periods.",
                         scm_fmt(o)), call. = FALSE)
    }
    invisible(TRUE)
}

# Parse a user-entered comma-separated list of numbers, e.g. "-2,-1,1,2".
# Duplicates are removed (they would otherwise produce duplicate table rows).
scm_parse_list <- function(s, name = "list") {
    parts <- trimws(strsplit(s, ",", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    v <- suppressWarnings(as.numeric(parts))
    if (length(v) == 0 || anyNA(v))
        stop(sprintf("'%s' must be a comma-separated list of numbers, e.g. -2,-1,1,2", name), call. = FALSE)
    unique(v)
}

# ---- panel construction and validation (used by the jamovi analysis and tests) --
# Converts a long data frame to the T x N outcome matrix and checks the design.
# Errors are signalled with stop(); the analysis converts them into jamovi messages.
scm_build_panel <- function(df, Y_name, unit_name, time_name, treated, treatmentTime,
                            exclude = character(0)) {
    if (anyNA(df[[unit_name]]))
        stop("The unit variable has missing values. Remove or fill in rows with a blank unit.", call. = FALSE)
    if (anyNA(df[[time_name]]))
        stop("The time variable has missing values. Remove or fill in rows with a blank time.", call. = FALSE)

    # jamovi may read the time variable as nominal (factor/text).
    df[[time_name]] <- suppressWarnings(as.numeric(as.character(df[[time_name]])))
    if (anyNA(df[[time_name]]))
        stop("The time variable must be numeric (e.g. years such as 1970, 1971, ...).", call. = FALSE)

    if (length(exclude) > 0) {
        if (treated %in% exclude)
            stop("The treated unit cannot appear in the donor exclusion list.", call. = FALSE)
        bad <- setdiff(exclude, unique(as.character(df[[unit_name]])))
        if (length(bad) > 0)
            stop(paste0("Unknown unit(s) in the exclusion list: ", paste(bad, collapse = ", ")), call. = FALSE)
        df <- df[!(as.character(df[[unit_name]]) %in% exclude), , drop = FALSE]
    }

    units <- unique(as.character(df[[unit_name]]))
    if (length(units) < 3)
        stop("Need at least two donor units after exclusions.", call. = FALSE)
    if (!(treated %in% units))
        stop("The treated unit was not found in the unit variable.", call. = FALSE)
    times <- sort(unique(df[[time_name]]))

    # Duplicate unit-period rows would be silently overwritten below.
    if (anyDuplicated(paste(df[[unit_name]], df[[time_name]])))
        stop("Each unit must have exactly one row per time period (duplicate unit-time rows found).", call. = FALSE)

    Y <- matrix(NA_real_, nrow = length(times), ncol = length(units),
                dimnames = list(as.character(times), units))
    Y[cbind(match(df[[time_name]], times),
            match(as.character(df[[unit_name]]), units))] <- df[[Y_name]]

    if (anyNA(Y))
        stop("Missing outcome values: the panel must be balanced, with an outcome for every unit in every period.", call. = FALSE)

    pre  <- which(times <  treatmentTime)
    post <- which(times >= treatmentTime)
    if (length(pre) < 2 || length(post) < 1)
        stop("Need at least two pre-treatment and one post-treatment period.", call. = FALSE)

    list(Y = Y, pre = pre, post = post, units = units, times = times, df = df)
}

# ---- slim states for the figures (only what each figure draws) -------------
scm_paths_state <- function(fit) {
    list(years = as.numeric(names(fit$synth)), treated = fit$treated_y,
         synth = fit$synth, post1 = fit$post[1])
}

# placebo: NULL or the slim object from scm_placebo_slim()
scm_gaps_state <- function(fit, placebo = NULL, k = 2) {
    list(years = as.numeric(names(fit$synth)), gap = fit$gap, post1 = fit$post[1],
         rmspe_pre = fit$rmspe_pre,
         gaps = if (is.null(placebo)) NULL else placebo$gaps,
         pre_rmspe = if (is.null(placebo)) NULL else placebo$placebo_pre_rmspe,
         k = k)
}

scm_weights_state <- function(fit) list(weights = fit$weights)
