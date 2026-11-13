# jamovi analysis class for the synthetic control method.
# Estimation lives in engine.R (scm_fit, scm_placebo, scm_loo, ...) so it stays
# testable without jamovi. This file wires options to results and states.

# Human-readable label for the p-value statistic.
meth_label <- function(m) if (identical(m, "att")) "two-sided, |average gap|" else "one-sided, RMSPE ratio"

# Turn engine errors into jamovi messages.
scm_guard <- function(expr) tryCatch(expr, error = function(e) jmvcore::reject(conditionMessage(e)))

scmClass <- R6::R6Class(
    "scmClass",
    inherit = scmBase,
    private = list(

        # ---- helpers -------------------------------------------------------
        .list = function(s) {
            v <- trimws(strsplit(if (is.null(s)) "" else s, ",", fixed = TRUE)[[1]])
            v[nzchar(v)]
        },

        # Donor names (after exclusions and without the treated unit), for .init.
        .donors = function(units, opts) {
            setdiff(units[!is.na(units)], c(opts$treatedUnit, private$.list(opts$donorExclude)))
        },

        # Write (value rows) into a robustness table keyed by its label column.
        .fillRows = function(tbl, labels, att, ratio, col) {
            tbl$deleteRows()
            for (i in seq_along(labels)) {
                vals <- list(att = att[i], ratio = ratio[i])
                vals[[col]] <- labels[i]
                tbl$addRow(rowKey = labels[i], values = vals)
            }
        },

        .fillMain = function(fit, opts, units, pre, post) {
            tbl <- self$results$summary
            tbl$deleteRows()
            tbl$addRow(rowKey = "pre",   values = list(metric = "RMSPE, pre-treatment",  value = fit$rmspe_pre))
            tbl$addRow(rowKey = "post",  values = list(metric = "RMSPE, post-treatment", value = fit$rmspe_post))
            tbl$addRow(rowKey = "ratio", values = list(metric = "Post / pre RMSPE ratio", value = fit$ratio))
            tbl$addRow(rowKey = "att",   values = list(metric = "Average post-treatment gap (ATT)",
                                                       value = scm_att(fit)))
            tbl$setNote("sample", sprintf(
                "Donor pool: %d units. Pre-treatment periods: %d. Post-treatment periods: %d.",
                length(units) - 1L, length(pre), length(post)))

            wt <- self$results$weights
            wt$deleteRows()
            w <- sort(fit$weights, decreasing = TRUE)
            for (u in names(w))
                wt$addRow(rowKey = u, values = list(unit = u, weight = w[[u]]))

            bal <- fit$balance
            bt <- self$results$balance
            bt$deleteRows()
            if (length(opts$predictors) == 0)
                bal$predictor <- paste(opts$outcome, bal$predictor)
            for (i in seq_len(nrow(bal)))
                bt$addRow(rowKey = bal$predictor[i], values = list(
                    predictor = bal$predictor[i], treated = bal$treated[i], synthetic = bal$synthetic[i]))
        },

        .fillPermutation = function(plc, method, k) {
            pt <- self$results$permTable
            pt$deleteRows()
            p_all <- scm_pvalue(plc, method = method)
            p_fil <- scm_pvalue(plc, max_pre_ratio = k, method = method)
            pt$addRow(rowKey = "all", values = list(sample = "All placebos",
                n = as.integer(p_all[["n_placebo"]]), p = p_all[["p"]]))
            pt$addRow(rowKey = "filtered", values = list(
                sample = sprintf("Pre-fit RMSPE ≤ %s × treated", scm_fmt(k)),
                n = as.integer(p_fil[["n_placebo"]]), p = p_fil[["p"]]))
            pt$setNote("method", sprintf(paste(
                "Statistic: %s. p is the share of placebos, counting the treated unit,",
                "with a statistic at least as large as the treated unit's."), meth_label(method)))
        },

        # Set up every row whose label is known before estimation.
        .init = function() {
            opts <- self$options
            tbl <- self$results$summary
            tbl$addRow(rowKey = "pre",   values = list(metric = "RMSPE, pre-treatment"))
            tbl$addRow(rowKey = "post",  values = list(metric = "RMSPE, post-treatment"))
            tbl$addRow(rowKey = "ratio", values = list(metric = "Post / pre RMSPE ratio"))
            tbl$addRow(rowKey = "att",   values = list(metric = "Average post-treatment gap (ATT)"))

            pt <- self$results$permTable
            pt$addRow(rowKey = "all", values = list(sample = "All placebos"))
            pt$addRow(rowKey = "filtered", values = list(sample = "Pre-fit filtered"))

            for (p in opts$predictors)
                self$results$balance$addRow(rowKey = p, values = list(predictor = p))

            t0 <- opts$treatmentTime
            add_scen <- function(tbl, text, labfun) {
                v <- tryCatch(scm_parse_list(text, "x"), error = function(e) NULL)
                if (is.null(v) || is.null(t0)) return(invisible(NULL))
                for (x in v) {
                    lab <- labfun(x)
                    tbl$addRow(rowKey = lab, values = list(scenario = lab))
                }
            }
            add_scen(self$results$timingTable,   opts$timingOffsets,   function(o) scm_lab_timing(t0 + o))
            add_scen(self$results$preLengthTable, opts$preLengthDrops, scm_lab_prelength)
            add_scen(self$results$intimeTable,   opts$intimeOffsets,   function(o) scm_lab_intime(t0 + o))

            if (!is.null(opts$treatedUnit) && !is.null(opts$unit) && !is.null(self$data)) {
                units <- tryCatch(unique(as.character(self$data[[opts$unit]])), error = function(e) NULL)
                if (!is.null(units)) {
                    donors <- private$.donors(units, opts)
                    sub <- private$.list(opts$looSubset)
                    if (length(sub) > 0) donors <- intersect(sub, donors)
                    for (u in donors)
                        self$results$looTable$addRow(rowKey = u, values = list(removed = u))
                }
            }
            invisible(NULL)
        },

        .run = function() {
            opts <- self$options
            if (isTRUE(opts$showHelp))
                self$results$helpText$setContent(scm_help_html())
            if (is.null(opts$outcome) || is.null(opts$unit) || is.null(opts$time) ||
                is.null(opts$treatedUnit)) {
                self$results$notes$setContent(
                    "<p>Select an outcome, a unit identifier, a time variable and a treated unit to begin.</p>")
                return()
            }
            self$results$notes$setContent("")

            nst <- as.integer(opts$optStarts)
            sd  <- as.integer(opts$seed)
            ck  <- function() private$.checkpoint()

            # ---- panel and inputs -------------------------------------------
            panel <- scm_guard(scm_build_panel(self$data, opts$outcome, opts$unit, opts$time,
                                               opts$treatedUnit, opts$treatmentTime,
                                               exclude = private$.list(opts$donorExclude)))
            df <- panel$df; Y <- panel$Y; pre <- panel$pre; post <- panel$post
            units <- panel$units; times <- panel$times
            t0 <- opts$treatmentTime

            # X_fun(pre_times) builds the predictor matrix from the given pre-treatment
            # periods only, so robustness windows never use their post-treatment data.
            X_fun <- function(pre_times)
                scm_predictor_matrix(df, opts$unit, opts$time, opts$predictors, units,
                                     pre_times, stat = opts$predictorStat)
            X <- scm_guard(X_fun(times[pre]))

            # Robustness settings are parsed and checked before anything is estimated.
            want <- function(a, b) isTRUE(opts[[a]]) || isTRUE(opts[[b]])
            timing <- if (want("timingCheck", "timingFig"))
                scm_guard(scm_parse_list(opts$timingOffsets, "Timing offsets")) else NULL
            drops <- if (want("preLengthCheck", "preLengthFig"))
                scm_guard(scm_parse_list(opts$preLengthDrops, "Pre-period drops")) else NULL
            intime <- if (want("intimeCheck", "intimeFig"))
                scm_guard(scm_parse_list(opts$intimeOffsets, "In-time placebo years")) else NULL
            scm_guard(scm_check_scenarios(times, t0, timing = timing, drops = drops, intime = intime))

            only <- NULL
            if (want("looCheck", "looFig")) {
                only <- private$.list(opts$looSubset)
                if (length(only) == 0) only <- NULL else {
                    if (opts$treatedUnit %in% only)
                        jmvcore::reject("The treated unit cannot be in the leave-one-out subset.")
                    bad <- setdiff(only, units)
                    if (length(bad) > 0)
                        jmvcore::reject(paste0("Unknown unit(s) in the leave-one-out subset: ",
                                               paste(bad, collapse = ", ")))
                }
            }

            # ---- main fit ---------------------------------------------------
            fit <- scm_guard(scm_fit(Y, treat = opts$treatedUnit, pre = pre, post = post,
                                     X = X, starts = nst, seed = sd))
            private$.fillMain(fit, opts, units, pre, post)
            self$results$pathsPlot$setState(scm_paths_state(fit))
            self$results$weightsPlot$setState(scm_weights_state(fit))
            ck()   # the main results are now on screen; jamovi may interrupt from here

            # ---- in-space placebo test --------------------------------------
            k <- opts$placeboPreRatio
            plc <- NULL
            if (isTRUE(opts$placebo)) {
                plc <- self$results$placeboState$state
                if (is.null(plc)) {
                    plc <- scm_guard(scm_placebo_slim(scm_placebo(
                        Y, treat = opts$treatedUnit, pre = pre, post = post, X = X,
                        fit_true = fit, checkpoint = ck, starts = nst, seed = sd)))
                    self$results$placeboState$setState(plc)
                }
                private$.fillPermutation(plc, opts$pValueMethod, k)
            }
            self$results$gapsPlot$setState(scm_gaps_state(fit, placebo = plc, k = k))

            # ---- robustness checks ------------------------------------------
            att_main <- scm_att(fit)

            if (want("looCheck", "looFig"))
                private$.runBlock("looState", "looTable", "looPlot", "removed", isTRUE(opts$looCheck),
                    function() scm_guard(scm_loo(Y, opts$treatedUnit, pre, post, X = X, only = only,
                                                 starts = nst, seed = sd, checkpoint = ck)),
                    function(st) st$removed,
                    function(st) list(loo = st, att = att_main))

            if (want("timingCheck", "timingFig"))
                private$.runBlock("timingState", "timingTable", "timingPlot", "scenario", isTRUE(opts$timingCheck),
                    function() scm_guard(scm_timing(Y, opts$treatedUnit, t0, offsets = timing,
                                                    X_fun = X_fun, starts = nst, seed = sd, checkpoint = ck)),
                    function(st) st$scenario,
                    function(st) list(df = st, att = att_main, xlab = "Treatment year"))

            if (want("preLengthCheck", "preLengthFig"))
                private$.runBlock("preLengthState", "preLengthTable", "preLengthPlot", "scenario",
                    isTRUE(opts$preLengthCheck),
                    function() scm_guard(scm_prelength(Y, opts$treatedUnit, t0, drop = drops,
                                                       X_fun = X_fun, starts = nst, seed = sd, checkpoint = ck)),
                    function(st) st$scenario,
                    function(st) list(df = st, att = att_main, xlab = "Last pre-treatment year used"))

            if (want("intimeCheck", "intimeFig"))
                private$.runBlock("intimeState", "intimeTable", "intimePlot", "scenario",
                    isTRUE(opts$intimeCheck),
                    function() scm_guard(scm_intime(Y, opts$treatedUnit, t0, offsets = intime,
                                                    X_fun = X_fun, starts = nst, seed = sd, checkpoint = ck)),
                    function(st) st$scenario,
                    function(st) list(df = st, att = att_main, xlab = "Fake treatment year"))
        },

        # One robustness check: reuse the stored computation if its options did not change,
        # otherwise compute it (with checkpoints), then fill its table and figure.
        .runBlock = function(state_item, table_item, fig_item, label_col, show_table,
                             compute, labels_of, fig_state) {
            st <- self$results[[state_item]]$state
            if (is.null(st)) {
                st <- compute()
                self$results[[state_item]]$setState(st)
            }
            if (show_table)
                private$.fillRows(self$results[[table_item]], labels_of(st), st$att, st$ratio, label_col)
            self$results[[fig_item]]$setState(fig_state(st))
            invisible(NULL)
        },

        .plotPaths = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_paths(image$state, ylab = self$options$outcome) + ggtheme)
            TRUE
        },

        .plotGaps = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_gaps(image$state) + ggtheme)
            TRUE
        },

        .plotWeights = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_weights(image$state$weights) + ggtheme)
            TRUE
        },

        .plotLoo = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_loo(image$state$loo, image$state$att) + ggtheme)
            TRUE
        },

        .plotTiming = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_scenario(image$state$df, image$state$att, image$state$xlab) + ggtheme)
            TRUE
        },

        .plotPreLength = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_scenario(image$state$df, image$state$att, image$state$xlab) + ggtheme)
            TRUE
        },

        .plotIntime = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(scm_plot_scenario(image$state$df, image$state$att, image$state$xlab) + ggtheme)
            TRUE
        },

        # Hidden result items that only hold reusable computations.
        .plotNone = function(image, ...) FALSE
    )
)
