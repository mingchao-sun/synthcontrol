# synthcontrol 1.0.0 (2026-11-12)

First public release.

- Panel construction and validation moved into `scm_build_panel()`, with stability tests for missing values, text time variables, unusual names, empty data, duplicate rows and invalid treatment settings.
- Method notes can be displayed in the results.
- Documentation and metadata finalised for the jamovi library.

# synthcontrol 0.2.0 (2026-11-12)

- In-space placebo test with two permutation statistics (post/pre RMSPE ratio, one-sided; absolute average gap, two-sided).
- Pre-fit filter for placebo units (configurable multiple of the treated unit's pre-treatment RMSPE).
- Robustness checks, each with table and/or figure: leave-one-out (optional donor subset), treatment timing offsets, pre-period length, in-time placebo.
- Predictors recomputed within each robustness window, so no post-treatment data enter the predictors of a placebo window.
- Donor exclusion; predictor summary (mean, median, last pre-treatment value); reproducible random starts (seed, number of starts).
- Method notes can be shown in the results (option "Show method notes in results").
- Time variable coerced to numeric; unbalanced panels and duplicate unit-time rows are rejected with a message.

# synthcontrol 0.1.1

- Initial jamovi module: donor weights, predictor weights, pre/post RMSPE, weights and balance tables, paths and gap plots.
