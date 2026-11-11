# plots.R ------------------------------------------------------------------
# Base-graphics plotting functions for synthetic control results.
# Each function draws on the current device so it can be used both inside
# jamovi (Image renderers) and in standalone scripts.

scm_plot_paths <- function(fit, ylab = "Outcome", main = "Treated vs synthetic") {
    yrs <- as.numeric(names(fit$synth))
    rng <- range(c(fit$treated_y, fit$synth))
    graphics::plot(yrs, fit$treated_y, type = "l", lwd = 2, col = "black",
                   ylim = rng, xlab = "Time", ylab = ylab, main = main, las = 1)
    graphics::lines(yrs, fit$synth, lwd = 2, lty = 2, col = "#1f77b4")
    graphics::abline(v = yrs[fit$post[1]] - 0.5, lty = 3, col = "grey40")
    graphics::legend("bottomleft", legend = c("Treated", "Synthetic"),
                     lty = c(1, 2), lwd = 2, col = c("black", "#1f77b4"),
                     bty = "n")
    invisible(NULL)
}

# max_pre_ratio: only show placebos whose pre-treatment RMSPE is at most
# max_pre_ratio times the treated unit's (Abadie et al., 2010 practice).
scm_plot_gaps <- function(fit, placebo = NULL, max_pre_ratio = 2,
                          ylab = "Gap (treated - synthetic)", main = "Gap in outcome") {
    yrs <- as.numeric(names(fit$synth))
    keep <- rep(TRUE, 0)
    G <- NULL
    if (!is.null(placebo)) {
        keep <- placebo$placebo_pre_rmspe <= max_pre_ratio * fit$rmspe_pre
        G <- placebo$gaps[, keep, drop = FALSE]
    }
    rng <- range(c(fit$gap, as.vector(G)), na.rm = TRUE)
    graphics::plot(NA, xlim = range(yrs), ylim = rng, xlab = "Time", ylab = ylab,
                   main = main, las = 1)
    if (!is.null(G)) {
        for (j in seq_len(ncol(G)))
            graphics::lines(yrs, G[, j], col = "grey75", lwd = 1)
    }
    graphics::lines(yrs, fit$gap, lwd = 2, col = "black")
    graphics::abline(h = 0, col = "grey40")
    graphics::abline(v = yrs[fit$post[1]] - 0.5, lty = 3, col = "grey40")
    if (!is.null(placebo))
        graphics::legend("bottomleft", legend = c("Treated", "Placebo (pre-fit ok)"),
                         lwd = c(2, 1), col = c("black", "grey75"), bty = "n")
    invisible(NULL)
}

scm_plot_weights <- function(fit, top = 15, main = "Donor weights") {
    w <- sort(fit$weights, decreasing = TRUE)
    w <- utils::head(w, top)
    graphics::par(mar = c(4, 8, 3, 1))
    graphics::barplot(rev(w), horiz = TRUE, las = 1, col = "#1f77b4", border = NA,
                      xlab = "Weight", main = main, xlim = c(0, max(w) * 1.1))
    invisible(NULL)
}

# Leave-one-out dot plot: ATT when each donor is removed, with the full-sample ATT.
scm_plot_loo <- function(loo, att_main, main = "Leave-one-out ATT") {
    d <- loo[order(loo$att), ]
    graphics::par(mar = c(4, 9, 3, 1))
    graphics::dotchart(d$att, labels = d$removed, cex = 0.7, pch = 19,
                       xlab = "Average post-treatment gap", main = main)
    graphics::abline(v = att_main, lty = 2, col = "#d62728")
    graphics::abline(v = 0, col = "grey60")
    invisible(NULL)
}

# Scenario plot for timing, pre-period length and in-time placebo checks.
scm_plot_scenario <- function(df, att_main, xlab, main) {
    graphics::plot(df$x, df$att, type = "b", pch = 19, col = "black",
                   ylim = range(c(df$att, att_main, 0)),
                   xlab = xlab, ylab = "Average post-treatment gap", main = main, las = 1)
    graphics::abline(h = att_main, lty = 2, col = "#d62728")
    graphics::abline(h = 0, col = "grey60")
    invisible(NULL)
}
