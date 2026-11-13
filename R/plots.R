# plots.R ------------------------------------------------------------------
# ggplot2 figures for synthetic control results. Each function takes the slim
# state built in engine.R (scm_paths_state(), scm_gaps_state(), ...) and returns
# a ggplot object. In jamovi the render function adds the user's theme
# (ggtheme) and prints the plot; titles come from the results heading.

scm_plot_paths <- function(st, ylab = "Outcome") {
    n <- length(st$years)
    df <- data.frame(year = rep(st$years, 2), value = c(st$treated, st$synth),
                     series = factor(rep(c("Treated", "Synthetic"), each = n),
                                     levels = c("Treated", "Synthetic")))
    ggplot2::ggplot(df, ggplot2::aes(x = year, y = value, linetype = series)) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = st$years[st$post1] - 0.5, linetype = "dotted", colour = "grey40") +
        ggplot2::labs(x = "Time", y = ylab, linetype = NULL) +
        ggplot2::theme(legend.position = "bottom")
}

# Gap plot. Placebo gaps are drawn only for placebos whose pre-treatment RMSPE
# is at most k times the treated unit's (Abadie et al., 2010 practice).
scm_plot_gaps <- function(st, ylab = "Gap (treated - synthetic)") {
    p <- ggplot2::ggplot()
    cap <- NULL
    if (!is.null(st$gaps)) {
        keep <- st$pre_rmspe <= st$k * st$rmspe_pre
        G <- st$gaps[, keep, drop = FALSE]
        if (ncol(G) > 0) {
            gl <- data.frame(year = rep(st$years, ncol(G)), gap = as.vector(G),
                             unit = rep(seq_len(ncol(G)), each = length(st$years)))
            p <- p + ggplot2::geom_line(data = gl,
                        ggplot2::aes(x = year, y = gap, group = unit),
                        colour = "grey70", linewidth = 0.4)
        }
        cap <- sprintf("Grey lines: placebo gaps with pre-fit RMSPE <= %g x treated (%d shown).",
                       st$k, ncol(G))
    }
    tr <- data.frame(year = st$years, gap = st$gap)
    p + ggplot2::geom_line(data = tr, ggplot2::aes(x = year, y = gap), linewidth = 0.9) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
        ggplot2::geom_vline(xintercept = st$years[st$post1] - 0.5, linetype = "dotted", colour = "grey40") +
        ggplot2::labs(x = "Time", y = ylab, caption = cap)
}

# Donor weights as horizontal bars (weights = named numeric vector).
scm_plot_weights <- function(weights, top = 15) {
    w <- utils::head(sort(weights, decreasing = TRUE), top)
    df <- data.frame(unit = factor(names(w), levels = rev(names(w))), weight = as.numeric(w))
    ggplot2::ggplot(df, ggplot2::aes(x = weight, y = unit)) +
        ggplot2::geom_col(fill = "#1f77b4") +
        ggplot2::labs(x = "Weight", y = NULL)
}

# Leave-one-out dot plot: ATT when each donor is removed, with the full-sample ATT.
scm_plot_loo <- function(loo, att_main) {
    d <- loo[order(loo$att), ]
    d$removed <- factor(d$removed, levels = d$removed)
    ggplot2::ggplot(d, ggplot2::aes(x = att, y = removed)) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
        ggplot2::geom_vline(xintercept = att_main, linetype = "dashed", colour = "#d62728") +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Average post-treatment gap", y = NULL,
                      caption = "Dashed line: full-sample ATT.")
}

# Scenario plot for timing, pre-period length and in-time placebo checks.
scm_plot_scenario <- function(df, att_main, xlab) {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = att)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
        ggplot2::geom_hline(yintercept = att_main, linetype = "dashed", colour = "#d62728") +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::labs(x = xlab, y = "Average post-treatment gap",
                      caption = "Dashed line: full-sample ATT.")
}
