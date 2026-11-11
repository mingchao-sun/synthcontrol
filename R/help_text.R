# Method notes shown in the results when 'Show method notes' is ticked.
# Generated from help/scm.md; regenerate this file if the help text changes.
scm_help_html <- function() r"---(
<h1>Synthetic Control</h1>
<p>Estimates the effect of a treatment on one unit (for example, a state that passed a policy) by comparing its outcome path with a weighted combination of untreated donor units. The method follows Abadie, Diamond and Hainmueller (2010).</p>
<h2>Data</h2>
<p>The data must be a <strong>balanced panel in long format</strong>: one row per unit and period. Every unit needs an outcome value in every period.</p>
<ul>
<li><strong>Outcome</strong>: numeric variable to analyse (for example, cigarette sales per capita).</li>
<li><strong>Unit</strong>: identifier of the cross-sectional unit (for example, state name).</li>
<li><strong>Time</strong>: period variable. Use numeric years (1970, 1971, ...) or consecutive period indices.</li>
<li><strong>Treated unit</strong>: the unit that receives the treatment. All other units form the donor pool.</li>
<li><strong>First treated period</strong>: the first period in which the treatment is in effect. Periods before this value are pre-treatment; periods from this value onward are post-treatment.</li>
</ul>
<h2>Predictors</h2>
<p>Optional. Each selected variable is summarised over the pre-treatment periods, separately for each unit (mean, median or last pre-treatment value; see <em>Predictor summary</em>). If no predictors are selected, the pre-treatment outcomes themselves are used as predictors. This is the outcome-only specification.</p>
<h2>Estimation</h2>
<p>Donor weights are non-negative and sum to one. For given predictor weights, they minimise the weighted discrepancy between the treated unit and the donor combination. The predictor weights are chosen to minimise the pre-treatment mean squared prediction error (MSPE) of the outcome. The optimisation over predictor weights is not globally convex, so it is restarted from several random starting points. Increase <em>Random starts</em> if results change across seeds.</p>
<ul>
<li><strong>Donor exclusion</strong>: drop units from the donor pool, for example states with confounding policies. The treated unit cannot be excluded.</li>
<li><strong>Random seed</strong>: makes the random starting points reproducible.</li>
</ul>
<h2>Results</h2>
<ul>
<li><strong>Summary</strong>: pre- and post-treatment RMSPE, their ratio, and the average post-treatment gap (ATT).</li>
<li><strong>Donor weights</strong>: weight on each donor, sorted from largest to smallest.</li>
<li><strong>Balance</strong>: treated unit and synthetic control values for each predictor (or pre-treatment outcome).</li>
<li><strong>Paths plot</strong>: treated and synthetic outcome over time.</li>
<li><strong>Gap plot</strong>: treated minus synthetic outcome over time. When the in-space placebo test is selected, grey lines show the placebo gaps that pass the pre-fit filter.</li>
<li><strong>Weights plot</strong>: donor weights as a bar chart.</li>
</ul>
<h2>In-space placebo test</h2>
<p>Each donor unit is treated as if it were the treated unit, and the synthetic control is re-estimated. The treated unit's effect is compared with this placebo distribution.</p>
<ul>
<li><strong>Permutation p-value statistic</strong>:<ul>
<li><em>RMSPE ratio (one-sided)</em>: ranks the post/pre RMSPE ratio. Appropriate when the hypothesis is a directional effect.</li>
<li><em>Absolute average gap (two-sided)</em>: ranks |ATT|. Appropriate when the direction of the effect is not specified in advance.</li>
</ul>
</li>
<li><strong>Pre-fit filter</strong>: placebo units whose pre-treatment RMSPE exceeds this multiple of the treated unit's pre-treatment RMSPE are excluded from the filtered p-value and from the gap plot. The default is 2. A poor pre-treatment fit makes the post/pre ratio uninformative.</li>
</ul>
<p>The summary table reports the p-value both with all placebos and with the filter applied.</p>
<h2>Robustness checks</h2>
<p>Each check can produce a table, a figure, or both.</p>
<ul>
<li><strong>Leave-one-out</strong>: re-estimates the model dropping each donor in turn. Optionally restrict to a comma-separated subset of donors. Large changes in the ATT indicate that a single donor drives the result.</li>
<li><strong>Treatment timing</strong>: re-estimates the model with the treatment date shifted by the offsets given (for example, -2, -1, 1, 2 years). Results should change little if the date is correctly chosen.</li>
<li><strong>Pre-period length</strong>: drops the given number of years from the end of the pre-treatment period. Shows sensitivity to the length of the pre-treatment window.</li>
<li><strong>In-time placebo</strong>: pretends the treatment occurred earlier, at the offsets given (for example, -6, -4, -2 years before treatment), using only pre-treatment data. A credible design should show effects that are small relative to the main ATT. Large placebo effects suggest that the pre-treatment fit does not capture the outcome dynamics well.</li>
</ul>
<h2>Notes and limitations</h2>
<ul>
<li>The estimator uses only the donor pool. Unmatched donors receive zero weight.</li>
<li>With few predictors relative to the number of donors, donor weights are generally not unique. Several weight vectors can match the treated unit equally well, and the solution then depends on the optimiser. Check the number of non-zero weights and the pre-treatment RMSPE before relying on the estimate, and prefer more informative predictors or the outcome-only specification.</li>
<li>Predictor-weight optimisation can take a minute or more with many random starts and many donors.</li>
<li>Permutation p-values have a minimum attainable value of 1 divided by the number of placebo units plus one. With few donors, the test has low power.</li>
<li>Results depend on the pre-treatment fit. Check the paths plot and the pre-treatment RMSPE before interpreting the post-treatment gap.</li>
</ul>
<h2>References</h2>
<p>Abadie, A., Diamond, A., and Hainmueller, J. (2010). Synthetic control methods for comparative case studies: estimating the effect of California's tobacco control program. <em>Journal of the American Statistical Association</em>, 105(490), 493–505.</p>
<p>Abadie, A., Diamond, A., and Hainmueller, J. (2011). Synth: An R package for synthetic control methods in comparative case studies. <em>Journal of Statistical Software</em>, 42(13), 1–17.</p>
<p>Abadie, A., Diamond, A., and Hainmueller, J. (2015). Comparative politics and the synthetic control method. <em>American Journal of Political Science</em>, 59(2), 495–510.</p>
<p>Abadie, A. (2021). Using synthetic controls: Feasibility, data requirements, and methodological aspects. <em>Journal of Economic Literature</em>, 59(2), 391–425.</p>
<p>Abadie, A., and Gardeazabal, J. (2003). The economic costs of conflict: A case study of the Basque Country. <em>American Economic Review</em>, 93(1), 113–132.</p>
<p>Ferman, B., and Pinto, C. (2021). Synthetic controls with imperfect pretreatment fit. <em>Quantitative Economics</em>, 12(4), 1197–1221.</p>
<p>Firpo, S., and Possebom, V. (2018). Synthetic control method: Inference, sensitivity analysis and confidence sets. <em>Journal of Causal Inference</em>, 6(2).</p>
)---"
