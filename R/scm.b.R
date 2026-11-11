
# jamovi analysis class for the synthetic control method.
# Heavy lifting is done in engine.R (scm_fit, scm_placebo, ...) so it stays testable.

# Human-readable label for the p-value method.
meth_label <- function(m) if (identical(m, "att")) "two-sided, |average gap|" else "one-sided, RMSPE ratio"

# Fill a results table from a robustness data frame (one row per scenario).
# label_name is the table's label column ("removed" for leave-one-out, else "scenario").
scm_fill_table <- function(tbl, df, scenario_col = "scenario", label_name = "scenario") {
    for (i in seq_len(nrow(df))) {
        key <- as.character(df[[scenario_col]][i])
        vals <- list(att = df$att[i], ratio = df$ratio[i])
        vals[[label_name]] <- key
        tbl$addRow(rowKey = key, values = vals)
    }
    invisible(NULL)
}

scmClass <- R6::R6Class(
    "scmClass",
    inherit = scmBase,
    private = list(

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

            Y_name <- opts$outcome
            unit_name <- opts$unit
            time_name <- opts$time
            excl <- trimws(strsplit(opts$donorExclude, ",", fixed = TRUE)[[1]])
            excl <- excl[nzchar(excl)]
            panel <- tryCatch(
                scm_build_panel(self$data, Y_name, unit_name, time_name, opts$treatedUnit,
                                opts$treatmentTime, exclude = excl),
                error = function(e) jmvcore::reject(conditionMessage(e)))
            df <- panel$df; Y <- panel$Y; pre <- panel$pre; post <- panel$post
            units <- panel$units; times <- panel$times

            # ---- predictors -------------------------------------------------
            # X_fun(pre_times) builds the K x N predictor matrix from the given
            # pre-treatment periods only. Robustness checks call it with their own
            # window, so predictors never use post-treatment data of that window.
            X_fun <- function(pre_times)
                scm_predictor_matrix(df, unit_name, time_name, opts$predictors, units,
                                     pre_times, stat = opts$predictorStat)
            X <- X_fun(times[pre])

            # ---- estimation settings ----------------------------------------
            nst <- as.integer(opts$optStarts)
            sd  <- as.integer(opts$seed)
            parse_or_reject <- function(s, name)
                tryCatch(scm_parse_list(s, name), error = function(e) jmvcore::reject(conditionMessage(e)))

            # ---- estimation -------------------------------------------------
            fit <- scm_fit(Y, treat = opts$treatedUnit, pre = pre, post = post, X = X,
                           starts = nst, seed = sd)

            # ---- summary table ----------------------------------------------
            tbl <- self$results$summary
            tbl$addRow(rowKey = "pre",   values = list(metric = "RMSPE, pre-treatment",  value = fit$rmspe_pre))
            tbl$addRow(rowKey = "post",  values = list(metric = "RMSPE, post-treatment", value = fit$rmspe_post))
            tbl$addRow(rowKey = "ratio", values = list(metric = "Post / pre RMSPE ratio", value = fit$ratio))
            tbl$addRow(rowKey = "att", values = list(metric = "Average post-treatment gap (ATT)",
                                                     value = scm_att(fit)))

            # ---- weights table ----------------------------------------------
            w <- sort(fit$weights, decreasing = TRUE)
            for (u in names(w))
                self$results$weights$addRow(rowKey = u, values = list(unit = u, weight = w[[u]]))

            # ---- balance table ----------------------------------------------
            bal <- fit$balance
            if (length(opts$predictors) == 0)
                bal$predictor <- paste(opts$outcome, bal$predictor)
            for (i in seq_len(nrow(bal)))
                self$results$balance$addRow(rowKey = i, values = list(
                    predictor = bal$predictor[i],
                    treated   = bal$treated[i],
                    synthetic = bal$synthetic[i]))

            # ---- figures and placebo inference ------------------------------
            self$results$pathsPlot$setState(fit)
            self$results$weightsPlot$setState(fit)
            k <- opts$placeboPreRatio
            if (isTRUE(opts$placebo)) {
                pl <- scm_placebo(Y, treat = opts$treatedUnit, pre = pre, post = post, X = X,
                                  fit_true = fit, starts = nst, seed = sd)
                meth <- opts$pValueMethod
                p_all <- scm_pvalue(pl, method = meth)
                p_fil <- scm_pvalue(pl, max_pre_ratio = k, method = meth)
                self$results$summary$addRow(rowKey = "p_all", values = list(
                    metric = sprintf("Permutation p-value, %s (all %d placebos)", meth_label(meth), as.integer(p_all[["n_placebo"]])),
                    value = p_all[["p"]]))
                self$results$summary$addRow(rowKey = "p_fil", values = list(
                    metric = sprintf("Permutation p-value, %s (pre-fit <= %g x treated, %d placebos)",
                                     meth_label(meth), k, as.integer(p_fil[["n_placebo"]])),
                    value = p_fil[["p"]]))
                self$results$gapsPlot$setState(list(fit = fit, placebo = pl, k = k))
            } else {
                self$results$gapsPlot$setState(list(fit = fit, placebo = NULL, k = k))
            }

            # ---- robustness checks (each: table and/or figure) ----------------
            att_main <- scm_att(fit)

            if (isTRUE(opts$looCheck) || isTRUE(opts$looFig)) {
                only <- trimws(strsplit(opts$looSubset, ",", fixed = TRUE)[[1]])
                only <- only[nzchar(only)]
                if (length(only) > 0) {
                    if (opts$treatedUnit %in% only)
                        jmvcore::reject("The treated unit cannot be in the leave-one-out subset.")
                    bad <- setdiff(only, units)
                    if (length(bad) > 0)
                        jmvcore::reject(paste0("Unknown unit(s) in the leave-one-out subset: ", paste(bad, collapse = ", ")))
                } else only <- NULL
                loo <- scm_loo(Y, opts$treatedUnit, pre, post, X = X, only = only, starts = nst, seed = sd)
                if (isTRUE(opts$looCheck))
                    scm_fill_table(self$results$looTable, loo, scenario_col = "removed", label_name = "removed")
                self$results$looPlot$setState(list(loo = loo, att = att_main))
            }

            if (isTRUE(opts$timingCheck) || isTRUE(opts$timingFig)) {
                tm <- scm_timing(Y, opts$treatedUnit, opts$treatmentTime,
                                 offsets = parse_or_reject(opts$timingOffsets, "Timing offsets"),
                                 X_fun = X_fun, starts = nst, seed = sd)
                if (isTRUE(opts$timingCheck))
                    scm_fill_table(self$results$timingTable, tm)
                self$results$timingPlot$setState(list(df = tm, att = att_main,
                    xlab = "Treatment year", main = "Treatment timing sensitivity"))
            }

            if (isTRUE(opts$preLengthCheck) || isTRUE(opts$preLengthFig)) {
                pl2 <- scm_prelength(Y, opts$treatedUnit, opts$treatmentTime,
                                     drop = parse_or_reject(opts$preLengthDrops, "Pre-period drops"),
                                     X_fun = X_fun, starts = nst, seed = sd)
                if (isTRUE(opts$preLengthCheck))
                    scm_fill_table(self$results$preLengthTable, pl2)
                self$results$preLengthPlot$setState(list(df = pl2, att = att_main,
                    xlab = "Last pre-treatment year used", main = "Pre-period length sensitivity"))
            }

            if (isTRUE(opts$intimeCheck) || isTRUE(opts$intimeFig)) {
                it2 <- scm_intime(Y, opts$treatedUnit, opts$treatmentTime,
                                  offsets = parse_or_reject(opts$intimeOffsets, "In-time placebo years"),
                                  X_fun = X_fun, starts = nst, seed = sd)
                if (isTRUE(opts$intimeCheck))
                    scm_fill_table(self$results$intimeTable, it2)
                self$results$intimePlot$setState(list(df = it2, att = att_main,
                    xlab = "Fake treatment year", main = "In-time placebo"))
            }

            self$results$notes$setContent(sprintf(
                "<p>Donor pool: %d units. Pre-treatment periods: %d. Post-treatment periods: %d.%s</p>",
                length(units) - 1L, length(pre), length(post),
                if (isTRUE(opts$placebo)) "" else " Tick 'In-space placebo test' to compute permutation p-values."))
        }
        ,

        .plotPaths = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            scm_plot_paths(image$state, ylab = self$options$outcome)
            TRUE
        },

        .plotGaps = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            scm_plot_gaps(image$state$fit, placebo = image$state$placebo,
                          max_pre_ratio = image$state$k)
            TRUE
        },

        .plotWeights = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            scm_plot_weights(image$state)
            TRUE
        }
        ,

        .plotLoo = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            scm_plot_loo(image$state$loo, image$state$att)
            TRUE
        },

        .plotTiming = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            st <- image$state
            scm_plot_scenario(st$df, st$att, st$xlab, st$main)
            TRUE
        },

        .plotPreLength = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            st <- image$state
            scm_plot_scenario(st$df, st$att, st$xlab, st$main)
            TRUE
        },

        .plotIntime = function(image, ...) {
            if (is.null(image$state)) return(FALSE)
            st <- image$state
            scm_plot_scenario(st$df, st$att, st$xlab, st$main)
            TRUE
        }
    )
)
