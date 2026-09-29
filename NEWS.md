# FastSurvival 1.0.0

## New features

* `coxph_fast()` gains a `strata` argument for a stratified Pike-Halley
  estimate of a common hazard ratio with a separate baseline hazard in each
  stratum (#1, suggested by Isaac Gravestock). The risk sets are formed within
  each stratum, and the observed and expected totals, the score, the
  information, and the curvature term are summed over strata, so the result
  approximates `coxph(... + strata(s), ties = "breslow")`. The `strata`
  argument of `analysis_fast()` now also stratifies `stat = "coxph"`.

## Bug fixes

* `simdata_fast()` with subgroups failed when the sample size was given as a
  scalar `n` with `alloc` (#2, reported by Isaac Gravestock), because a scalar
  `n` with a common prevalence was always treated as a single group. A scalar
  `n` now gives a two-group simulation when `alloc` is supplied explicitly or
  when a survival or dropout specification is a per-group list with per-cell
  elements (such as `list(list(0.10, 0.08, 0.06), 0.05)`), and the result is
  identical to the per-group `n` specification. Without either signal, a list
  with subgroups is still read as per-cell values of a single group. The rules
  are described in the Details of `?simdata_fast`. In the same area:
  * a shared (non-list) survival or dropout specification is accepted in a
    two-group simulation with subgroups;
  * lists of the wrong length are reported instead of being silently
    truncated or failing with "subscript out of bounds";
  * per-cell `e.time` and `d.time` lists are honored in a one-group
    simulation with subgroups;
  * group-specific prevalence must use the same number of levels for each
    factor in both groups;
  * the illness-death model also treats an explicit `alloc` as a two-group
    request;
  * without subgroups, a scalar `n` with an explicit `alloc` and a shared
    (non-list) hazard now gives two groups sharing that hazard, where it
    previously gave one group of size `n` and ignored `alloc`.
* `analysis_fast()` with `stat = "milestone"`, `ms.method = "loglog"`, and
  `side = 1` returned the upper-tail p-value although the log-log statistic
  is negative under treatment benefit, so the one-sided p-value was
  approximately one minus the correct value. It now uses the lower tail, as
  `milestone_fast()` does.
* `analysis_fast()` now reports the `"ahsw"` p-values according to `side`, as
  `ahsw_fast()` does; they were always two-sided before.
* `pairwise_fast()` in event-driven mode reported every administratively
  censored subject as a dropout (`n.pipeline` was always 0). Dropout and
  pipeline counts are now computed at the per-simulation cutoff. The
  event-driven mode also keeps the other columns of `data`, so `strata` works
  there, and an ambiguous `p.col` (several statistics) is reported clearly.
* The log-rank family (`survdiff_fast()`, `maxcombo_fast()`, `rmw_fast()`, and
  the weighted and stratified variants) left out of the observed and expected
  counts an event at a time when only one subject was at risk. The test
  statistics were unaffected, but the printed counts differed from
  `survival::survdiff()`.
* `medsurv_fast()` (and `stat = "medsurv"` in `analysis_fast()`) now follows
  the `survival::survfit()` median convention: the comparison with 0.5 uses a
  tolerance, and a curve that equals 0.5 on a flat stretch gives the midpoint
  of that stretch. Before, floating-point rounding could move the median to
  the next event time. With `method = "km"` the kernel hazard is evaluated at
  the new median, so its standard error changes in the flat-stretch case as
  well. The printed median of `kmcurve_fast()` objects uses the
  same step-function rule instead of linear interpolation.
* `simdata_fast()` with a scalar `n` split the total with `round()`, which
  could lose or add a subject (for example `n = 7` gave 4 + 4). The split now
  always adds up to `n`, and `alloc` is validated.
* `simdata_fast()` with `fixed.alloc = TRUE` assigned the subgroup cells in
  contiguous blocks, so subgroup membership was tied to the accrual interval.
  The fixed labels are now randomly permuted within each simulation, which
  changes the generated data for `fixed.alloc = TRUE` (only) for a given seed.
  The fixed counts are also protected against floating-point shares such as
  `100 * 0.29`.
* `simdata_fast()` in the illness-death model now treats a per-group
  `d.hazard` or `d.median` list as a two-group request, as in the
  single-endpoint model.
* `milestone_fast()` (and `stat = "milestone"`) returned a `NaN` standard
  error when a Kaplan-Meier curve reached zero by the milestone; it is now 0,
  as in `survfit_fast()`.
* `survfit_fast()` caps the upper limit of the `"log"` interval at 1 and gives
  a degenerate interval when the standard error is zero, as
  `survival::survfit()` does.
* The print method of `survdiff_fast()` labels the unweighted stratified
  test as stratified.
* `simsummary_fast()` silently kept only the last row when `data` had more
  than one row per simulation and look, such as the stacked arms of
  `pairwise_fast()` output. It now summarizes each arm separately when `data`
  has an `arm` column (the output gains an `arm` column and the print method a
  heading per arm), and stops with an error for other duplicated rows.
* The print method of `simsummary_fast()` no longer runs the label "Expected
  analysis time at stop:" into its value.

## Input validation

* `rmst_fast()`, `ahsw_fast()`, `milestone_fast()`, `wmst_fast()` (a supplied
  `tau2`), and `ahr_fast()` (a supplied `tau`) warn when the truncation time
  or milestone exceeds the largest observed time of a group, where the
  Kaplan-Meier curve is not estimated and is carried forward flat.
* `survdiff_fast()`, `coxph_fast()`, `rmst_fast()`, `maxcombo_fast()`,
  `rmw_fast()`, `ahsw_fast()`, `survfit_fast()`, and `analysis_fast()` now
  check that `event` is coded 0/1 without missing values and, for the
  two-group functions, that `group` has exactly two values without missing
  values and that `control` is one of them. Before, a wrong `control` label or a
  1/2 event coding silently produced wrong or missing results.
* `analysis_fast()` requires whole-number `event.looks` and recodes factor and
  character subgroup columns consistently, so `by.subgroup = TRUE` labels the
  populations correctly for such columns.
* `survfit_fast()` checks that `t_sorted` and `e_sorted` have the same length
  and, with `presorted = TRUE`, that the times are sorted.
* `simdata_fast()` checks that hazards are non-negative and that piecewise
  breakpoints start at 0 and increase, and that the output fits in an R
  vector.
* `milestone_fast()`, `medsurv_fast()`, and `wmst_fast()` reject missing
  times or groups.

## Documentation

* The modestly-weighted weight cap is documented as `1 / S(t_star-)`, the
  pooled Kaplan-Meier value just before `t_star`, which is what the code
  computes (as in nphRCT).
* The speed-comparison vignette reports speed gains re-measured for this
  release.
* The validation vignette compares the stratified `coxph_fast()` with
  `coxph()` using `strata()`.
* The `ahr_fast()` documentation explains that, as in the `AHR` package, the
  estimate is not symmetric in the groups when both groups have events at the
  same time, so `control` should be the actual reference group.
* The `ahsw_fast()` documentation states the actual condition for `NA`
  results (no events up to `tau` in a group).
* The `analysis_fast()` documentation no longer lists a `look.type` column,
  which the function does not return.

# FastSurvival 0.2.0

* New estimation and testing functions:
  * `rmst_fast()`: restricted mean survival time for a single group or a
    two-group comparison (difference and ratio contrasts), integrating the
    Kaplan-Meier survival step function in a single C++ scan.
  * `wmst_fast()`: window mean survival time over an interval, generalizing
    `rmst_fast()` (which is the special case with a lower window limit of
    zero), for a single group or a two-group difference, computed in the same
    single C++ scan with a Greenwood-type variance in which each event time
    contributes its squared remaining window area.
  * `milestone_fast()`: two-group comparison of Kaplan-Meier survival at a
    milestone timepoint, with Wald, log-log, and MOVER inference methods.
  * `medsurv_fast()`: median survival time for a single group or a two-group
    difference, with a native kernel-hazard variance method and an
    `nph`-compatible local-constant-hazard method that reproduces the median
    comparison of the `nph` package to numerical precision; the point estimate
    is the same under both methods.
  * `maxcombo_fast()`: max-combo test over a set of Fleming-Harrington
    weighted log-rank statistics, with the joint p-value obtained from the
    implied multivariate normal distribution.
  * `rmw_fast()`: robust modestly-weighted log-rank test of Magirr and Öhrn,
    the maximum of the standard log-rank and a modestly-weighted log-rank
    statistic, with the joint p-value obtained from the implied bivariate
    normal distribution.
  * `wkm_fast()`: weighted Kaplan-Meier (Pepe-Fleming) test, the weighted
    integrated difference between two Kaplan-Meier curves, with Pepe-Fleming,
    square-root, and constant weights, reproducing the weighted Kaplan-Meier
    statistic of the `nphsim` package.
  * `ahsw_fast()`: average hazard with survival weight of Uno and Horiguchi,
    reporting the ratio (RAH) and difference (DAH) contrasts.
  * `ahr_fast()`: Kalbfleisch-Prentice average hazard ratio between two groups
    over a restricted interval, the estimator used by Dormuth et al. (2024)
    for sample-size calculation under non-proportional hazards, with a test on
    the group-share scale and an equivalent test and confidence interval on
    the log scale.
* `survdiff_fast()` gains weighted log-rank tests (Fleming-Harrington,
  modestly-weighted, Gehan-Breslow, Tarone-Ware) and stratified and
  stratified-weighted variants, all sharing the single-scan C++ backend.
* New simulation layer:
  * `simdata_fast()` extended with optional subgroups defined by a prevalence
    specification and a flexible accrual specification: `a.rate` gives absolute
    accrual rates (with the end of an open final interval solved from the total
    when a trailing rate is supplied) and `a.prop` gives accrual proportions,
    with deterministic per-interval accrual counts. The entire generation
    pipeline runs in a single C++ kernel that materializes the output data
    frame once. It can also generate two correlated time-to-event
    endpoints (for example progression-free and overall survival) from an
    illness-death model with three transition hazards and optional
    treatment switching, reducing to the Fleischer maximal-independence
    model when the post-event hazard equals the direct terminal hazard.
    A vector `n` of length greater than two together with a per-arm
    survival list generates a multi-arm trial, each arm produced with the
    single-group kernel over a common accrual window and labeled 1 to
    `length(n)`, for analysis as pairwise contrasts against a shared control.
  * `analysis_fast()`: interim or sequential analysis of simulated data at one
    or more looks, defined by target event counts or calendar times, computed
    by a fused C++ kernel that reuses the analysis cores of the standalone
    functions. Supports subgroup analyses.
  * `pairwise_fast()`: runs `analysis_fast()` for each experimental arm against a
    shared control on multi-arm data, at either fixed calendar looks or the
    per-simulation cutoffs of a designated primary contrast, and stacks the
    results with an optional Bonferroni adjustment across contrasts.
  * `simsummary_fast()`: operating-characteristic summary (rejection and
    futility rates, stopping-look distribution, expected timing) from
    `analysis_fast()` output and supplied group-sequential boundaries, with a
    `print()` method that lays the results out as a group-sequential design
    report.
* New visualization layer:
  * `gen_scenario_fast()`: assembles one or more two-group scenarios into a
    `scenario_fast` object for design-stage exploration, with a `plot()`
    method that draws the analytic survival curves and the piecewise hazard
    ratio of each scenario and a `print()` method that summarizes the
    medians, the start and end hazard ratios, and whether the curves cross.
  * `kmcurve_fast()`: builds the Kaplan-Meier curves of a single trial
    realization (for example one replicate of `simdata_fast()`) into a
    `kmcurve_fast` object, with a `plot()` method that draws the curves with
    optional restricted-mean shading and a smoothed time-varying hazard-ratio
    panel, and a `print()` method that summarizes the events and medians.
* Each estimation and testing function has a corresponding `print()` method,
  and the print methods share a unified display format.
* New vignettes accompany the analysis, simulation, and visualization layers:
  validation against established packages, a speed comparison, a group-sequential
  design reproduction, a log-rank versus RMST comparison under nonproportional
  hazards, the Freidlin-Korn strong-null investigation, a correlated PFS and OS
  group-sequential design under the Fleischer model, a multiregional
  regional-consistency evaluation, and a multi-arm design analyzed as pairwise
  contrasts against a shared control.

# FastSurvival 0.1.0

* Initial release.
* Core computations implemented in `C++` via `Rcpp` for use inside large
  simulation loops.
* `survfit_fast()`: single-time-point Kaplan-Meier estimator with Greenwood
  standard error and plain / log / log-log confidence intervals. The C++
  backend locates the evaluation cutoff via binary search and accumulates
  the Kaplan-Meier product and Greenwood variance sum in a single scan over
  event positions. Returns an object of class `"survfit_fast"` with a
  `print()` method.
* `survdiff_fast()`: log-rank test returning a one-sided Z-score or a
  two-sided chi-square statistic. The C++ backend uses a two-pointer merge
  scan over pooled sorted vectors, eliminating the rank construction,
  `tabulate()`, and reverse cumulative sum operations of the standard
  implementation. Returns an object of class `"survdiff_fast"` with a
  `print()` method.
* `coxph_fast()`: closed-form hazard ratio estimator via the Pike-Halley
  Estimator method with Wald confidence interval. The C++ backend performs
  group splitting, at-risk counting, and per-distinct-event-time
  accumulation in a single pass. Returns an object of class `"coxph_fast"`
  with a `print()` method.
* `simdata_fast()`: clinical trial data simulator supporting one- and
  two-group designs, piecewise uniform accrual, and simple and piecewise
  exponential survival and dropout times. C++ backends handle piecewise
  sampling and two-group interleaving, and random number generation uses
  `dqrng`.
