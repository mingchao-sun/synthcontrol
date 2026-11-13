# Synthetic Control

Estimates the effect of a treatment on one unit (for example, a state that passed a policy) by comparing its outcome path with a weighted combination of untreated donor units. The method follows Abadie, Diamond and Hainmueller (2010). The references are listed under *References* in the results.

## Data

The data must be a **balanced panel in long format**: one row per unit and period. Every unit needs an outcome value in every period.

- **Outcome**: numeric variable to analyse (for example, cigarette sales per capita).
- **Unit**: identifier of the cross-sectional unit (for example, state name).
- **Time**: period variable. Use numeric years (1970, 1971, ...) or consecutive period indices.
- **Treated unit**: the unit that receives the treatment. All other units form the donor pool.
- **First treated period**: the first period in which the treatment is in effect. Periods before this value are pre-treatment; periods from this value onward are post-treatment.

## Predictors

Optional. Each selected variable is summarised over the pre-treatment periods, separately for each unit (mean, median or last pre-treatment value; see *Predictor summary*). If no predictors are selected, the pre-treatment outcomes themselves are used as predictors. This is the outcome-only specification.

## Estimation

Donor weights are non-negative and sum to one. For given predictor weights, they minimise the weighted discrepancy between the treated unit and the donor combination. The predictor weights are chosen to minimise the pre-treatment mean squared prediction error (MSPE) of the outcome. The optimisation over predictor weights is not globally convex, so it is restarted from several random starting points. Increase *Random starts* if results change across seeds.

- **Donor exclusion**: drop units from the donor pool, for example states with confounding policies. The treated unit cannot be excluded.
- **Random seed**: makes the random starting points reproducible.

## Results

- **Summary**: pre- and post-treatment RMSPE, their ratio, and the average post-treatment gap (ATT).
- **Donor weights**: weight on each donor, sorted from largest to smallest.
- **Balance**: treated unit and synthetic control values for each predictor (or pre-treatment outcome).
- **Paths plot**: treated and synthetic outcome over time.
- **Gap plot**: treated minus synthetic outcome over time. When the in-space placebo test is selected, grey lines show the placebo gaps that pass the pre-fit filter.
- **Weights plot**: donor weights as a bar chart.

## In-space placebo test

Each donor unit is treated as if it were the treated unit, and the synthetic control is re-estimated. The treated unit's effect is compared with this placebo distribution.

- **Permutation p-value statistic**:
    - *RMSPE ratio (one-sided)*: ranks the post/pre RMSPE ratio. Appropriate when the hypothesis is a directional effect.
    - *Absolute average gap (two-sided)*: ranks |ATT|. Appropriate when the direction of the effect is not specified in advance.
- **Pre-fit filter**: placebo units whose pre-treatment RMSPE exceeds this multiple of the treated unit's pre-treatment RMSPE are excluded from the filtered p-value and from the gap plot. The default is 2. A poor pre-treatment fit makes the post/pre ratio uninformative.

The summary table reports the p-value both with all placebos and with the filter applied.

## Robustness checks

Each check can produce a table, a figure, or both.

- **Leave-one-out**: re-estimates the model dropping each donor in turn. Optionally restrict to a comma-separated subset of donors. Large changes in the ATT indicate that a single donor drives the result.
- **Treatment timing**: re-estimates the model with the treatment date shifted by the offsets given (for example, -2, -1, 1, 2 years). Results should change little if the date is correctly chosen.
- **Pre-period length**: drops the given number of years from the end of the pre-treatment period. Shows sensitivity to the length of the pre-treatment window.
- **In-time placebo**: pretends the treatment occurred earlier, at the offsets given (for example, -6, -4, -2 years before treatment), using only pre-treatment data. A credible design should show effects that are small relative to the main ATT. Large placebo effects suggest that the pre-treatment fit does not capture the outcome dynamics well.

## Notes and limitations

- The estimator uses only the donor pool. Unmatched donors receive zero weight.
- With few predictors relative to the number of donors, donor weights are generally not unique. Several weight vectors can match the treated unit equally well, and the solution then depends on the optimiser. Check the number of non-zero weights and the pre-treatment RMSPE before relying on the estimate, and prefer more informative predictors or the outcome-only specification.
- Predictor-weight optimisation can take a minute or more with many random starts and many donors.
- Permutation p-values have a minimum attainable value of 1 divided by the number of placebo units plus one. With few donors, the test has low power.
- Results depend on the pre-treatment fit. Check the paths plot and the pre-treatment RMSPE before interpreting the post-treatment gap.
