# synthcontrol: synthetic control method for jamovi

Synthetic control method (Abadie, Diamond and Hainmueller, 2010) for comparative
case studies with panel data, as a jamovi module.

## Features

- Donor weights (non-negative, sum to one) and predictor weights, with pre- and post-treatment RMSPE.
- Optional predictors (mean, median or last pre-treatment value); outcome-only specification if none are given.
- In-space placebo test with one-sided (RMSPE ratio) or two-sided (absolute average gap) permutation p-values, with a pre-fit filter.
- Robustness checks: leave-one-out, treatment timing, pre-period length, in-time placebo. Each produces a table, a figure, or both.
- Donor exclusion and reproducible optimisation (random seed and starts).

## Installation

Tested with jamovi 28.17 (jmvtools 28.17). Requirements: R with the
packages `jmvcore`, `R6` and `quadprog`.

```r
# from a source folder, with jmvtools installed
jmvtools::install(pkg = "path/to/synthcontrol", home = "path/to/jamovi")
```

## Data

A balanced panel in long format: one row per unit and period, with numeric
time (years or consecutive periods) and a numeric outcome. The example data
`prop99_long.csv` (cigarette sales, 39 US states, 1970–2000) is included.

## Tests

Run from the module folder (the tests read `data/` by relative path):

```r
Rscript tests/run_tests.R
Rscript tests/test_robustness.R
Rscript tests/test_stability.R
```

## Help and limitations

Method notes are in `help/scm.md`, and can be displayed in the results with
the option *Show method notes in results*.

Donor weights are not unique when there are few predictors relative to the
number of donors. Check the number of non-zero weights and the pre-treatment
RMSPE before interpreting the estimate. Permutation p-values have a minimum of
1/(n+1) and low power with few donors.

## Citation

Abadie, A., Diamond, A., and Hainmueller, J. (2010). Synthetic control methods
for comparative case studies: estimating the effect of California's tobacco
control program. *Journal of the American Statistical Association*, 105(490), 493–505.

## Licence

GPL-3. See `LICENSE`.
