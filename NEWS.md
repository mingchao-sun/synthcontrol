# synthcontrol 1.0.0 (2026-11-12)

First public release.

- Panel construction and validation moved into `scm_build_panel()`, with stability tests for missing values, text time variables, unusual names, empty data, duplicate rows and invalid treatment settings.
- Method notes can be displayed in the results (option "Method notes", in the Output box).
- Long computations (placebo test and robustness checks) yield to jamovi between fits, and the placebo test and each robustness check are stored and reused until an option they depend on changes.
- Missing values in the outcome, predictors, unit or time variables give specific messages naming the variable, predictor and units.
- Robustness settings (timing offsets, pre-period drops, in-time placebo offsets, leave-one-out subset) are checked against the data before estimation, with messages naming the offending value.
- Permutation p-values are shown in their own table, with the statistic in a footnote. Table rows are set up before the analysis runs.
- Results carry jamovi references (Abadie et al., 2010; quadprog). Result tables and figures use jamovi's clearWith.
- Options panel: Placebo Test and Robustness Checks start collapsed; dependent options are disabled until their check is ticked; count and seed options are whole numbers.
- The Prop 99 example dataset is listed in jamovi's Data Library.
- Figures are ggplot2 plots and take the user's jamovi plot theme.
- The test suite is in testthat (`tests/testthat/`, run with `Rscript tests/run_all.R`).
- Documentation and metadata finalised for the jamovi library.

# synthcontrol 0.2.0 (2026-11-12)

- In-space placebo test with two permutation statistics (post/pre RMSPE ratio, one-sided; absolute average gap, two-sided).
- Pre-fit filter for placebo units (configurable multiple of the treated unit's pre-treatment RMSPE).
- Robustness checks, each with table and/or figure: leave-one-out (optional donor subset), treatment timing offsets, pre-period length, in-time placebo.
- Predictors recomputed within each robustness window, so no post-treatment data enter the predictors of a placebo window.
- Donor exclusion; predictor summary (mean, median, last pre-treatment value); reproducible random starts (seed, number of starts).
- Method notes can be shown in the results (option "Show method notes in results"; renamed "Method notes" in 1.0.0).
- Time variable coerced to numeric; unbalanced panels and duplicate unit-time rows are rejected with a message.

# synthcontrol 0.1.1

- Initial jamovi module: donor weights, predictor weights, pre/post RMSPE, weights and balance tables, paths and gap plots.
