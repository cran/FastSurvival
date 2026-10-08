# FastSurvival 1.1.0

## New features

* New `cutoff_fast()` computes the calendar time of each analysis look in
  every simulated trial in a single C++ pass, from a target number of events, a
  planned calendar time, a maximum calendar time, a minimum time after the
  previous look, and a minimum follow-up after a given number of enrolled
  subjects, combined as `min(max(...), maximum time)` as in
  `simtrial::get_analysis_date()`. Events can be counted on a subset of the
  subjects (for example the control group) or on another endpoint (for example
  the PFS columns of an illness-death simulation), and the event targets can be
  given per simulated trial.

* `analysis_fast()` gains a `cutoff.looks` argument, a matrix of per-simulation
  calendar cutoffs (such as the output of `cutoff_fast()`), as an alternative
  to `event.looks` and `time.looks`. This analyzes an endpoint at cutoffs
  determined by another endpoint, by a subset of the subjects, or by combined
  trigger rules, within the same fused C++ loop.

* `pairwise_fast()` gains a `cutoff.looks` argument, and its event-driven mode
  now computes the shared cutoffs with `cutoff_fast()` and analyzes every
  contrast through `cutoff.looks` instead of re-cutting the data in R. The
  results are unchanged.

* New `switch_fast()` applies treatment switching to simulated data, changing
  only the outcomes after each subject's switch: at the intermediate event of
  an illness-death simulation (for example progression), at an opening time
  after an interim analysis (crossover at a milestone, optionally only in the
  simulated trials selected by an interim decision), or at the later of the
  two. The remaining time to the terminal event is multiplied by an
  acceleration factor (the causal model of the RPSFT method) or redrawn from a
  new (piecewise) exponential hazard. Outcomes observed before the opening
  time are unchanged, so an interim analysis is unaffected by the switching.

* `analysis_fast()` gains an `mc.alpha` argument, the nominal levels at which
  the `"maxcombo"` p-values are compared. The p-value is then integrated only
  when the Bonferroni bounds `p_min <= p <= K * p_min` do not decide the
  comparison with the level; otherwise the bound on the same side of the level
  is reported. The decisions at the levels (for example in
  `simsummary_fast(p.col = "maxcombo.p", alpha = mc.alpha)`) are those of the
  exact p-values, and a `maxcombo.p.exact` column marks the integrated rows.
  The default (`NULL`) computes every p-value as before.

* The illness-death model of `simdata_fast()` now supports subgroups through
  `prevalence` and `fixed.alloc` (#3, suggested by Isaac Gravestock). The
  transition hazards, `switch.prop`, and the dropout specification accept
  per-group and per-cell nested lists with the same rules as `e.hazard`, for
  example `h01.hazard = list(list(0.10, 0.06), 0.05)`, so a mixture of
  illness-death models is simulated in one call. All subjects share one
  accrual process and the subgroups are assigned after accrual, so they enroll
  over the same calendar in proportion to their prevalence. Without
  `prevalence` the results are unchanged.

* `simdata_fast()` gains a `stream` argument that selects an independent
  `dqrng` random-number stream for the given `seed`, so that a large
  simulation can be generated in reproducible batches, sequentially or in
  parallel. Without `stream` the results are unchanged.

## Bug fixes

* `analysis_fast()` with `stat = "maxcombo"`, `side = 2`, and two or three
  weights failed with "TVPACK either needs all(lower == -Inf) or all(upper ==
  Inf)", because the two-sided rectangle was passed to the TVPACK algorithm,
  which handles only half-spaces. It now uses the GenzBretz algorithm in that
  case, as `maxcombo_fast()` does.

* `wkm_fast()` and `stat = "wkm"` of `analysis_fast()` divided the variance
  terms by the selected weight instead of the Pepe-Fleming weight, which is
  the censoring factor of the variance whatever weight is used. The results
  with `weight = "PF"` (the default) are unchanged. With `"sqrtPF"` and
  `"constant"` the standard error was too small and the test anti-conservative;
  for the constant weight under uniform censoring, the empirical standard
  deviation of the weighted difference was about 1.15 times the mean standard
  error. The error dates from version 0.2.0.

* `simdata_fast()` recorded a subject with neither a finite survival time nor
  a finite dropout time (possible with a zero hazard in the last piece, such
  as a cure fraction, and no dropout) as an event at an infinite time, so an
  event-driven look of `analysis_fast()` could be reached at an infinite
  cutoff. Such a subject now has `tte = Inf` and `event = 0`, in the
  single-endpoint and illness-death models and in `switch_fast()`, and is
  counted as in follow-up, not as a dropout, at an unreached look of
  `analysis_fast()`. A zero hazard no longer gives `NaN` when the draw equals
  the cumulative hazard at the last breakpoint.

* In the illness-death model of `simdata_fast()`, a switch through
  `switch.prop` could occur at an intermediate event after dropout, which is
  not observed. A switch now requires the intermediate event on or before
  dropout. The observed columns are unchanged; only the latent terminal time,
  `switched`, and `switch_time` of the affected subjects (who have dropped out
  before the intermediate event) change. The column `intermediate` describes
  the latent process, as the documentation now states.

* `milestone_fast()` and `stat = "milestone"` with `method = "loglog"` or
  `"mover"` returned missing confidence limits when the Kaplan-Meier estimate
  of a group was 0 or 1 at `tau`. The one-sample interval of that group now
  degenerates to the estimate, so the interval of the difference is reported;
  the `"loglog"` test statistic remains `NA` in that case.

* `simdata_fast()` with a single-level prevalence (for example
  `prevalence = 1`) did not output the subgroup column.

* `analysis_fast()` converted an `event.looks` target beyond the integer range
  to an undefined integer in C++. Such a target is now capped, as in
  `cutoff_fast()`, and is never reached.

* `switch_fast()` matched a named `cutoff` vector by position. The names are
  now matched to `data$sim`, as in the `cutoff.looks` argument of
  `analysis_fast()`.

## Input validation

* With `presorted = TRUE`, `survdiff_fast()`, `coxph_fast()`, `rmst_fast()`,
  `wmst_fast()`, `milestone_fast()`, `medsurv_fast()`, `maxcombo_fast()`,
  `rmw_fast()`, `wkm_fast()`, `ahsw_fast()`, and `ahr_fast()` now check the
  order in a single pass (with strata, also that the rows of each stratum are
  contiguous) and give an error for unsorted input instead of a wrong result,
  as `survfit_fast()` already did. The fused kernel of `analysis_fast()` is not
  affected.

* `simdata_fast()` checks that `n` contains positive whole numbers.

* `milestone_fast()`, `ahr_fast()`, and `kmcurve_fast()` check the event
  indicator before converting it to integer, so a value such as 0.7 is an
  error rather than a censoring, and `kmcurve_fast()` also rejects missing
  values. Negative times are rejected by `milestone_fast()`,
  `kmcurve_fast()`, and the functions that share the internal time and event
  check.

* `cutoff_fast()` warns when the cutoff of a look precedes that of the
  previous look in some simulated trials, for example when only the later look
  has a calendar cap.

## Documentation

* New vignette "Treatment switching and crossover after an interim analysis"
  simulates crossover after a positive PFS analysis and switching at
  progression, and checks the switching model with the RPSFT estimator of the
  rpsftm package, which is added to Suggests and used only when it is
  installed.

* The correlated PFS and OS vignette now analyzes OS at the PFS-driven cutoffs
  with `cutoff_fast()` and `cutoff.looks`.

* The validation vignette compares `cutoff_fast()` with
  `simtrial::get_analysis_date()` and the survival times generated by
  `switch_fast()` with their analytic distribution. It no longer states a
  reproduction of published p-values that the package sources do not contain,
  and it describes the p-values of `nph::logrank.maxtest()` correctly (they are
  two-sided).

* New vignette "Using your own data generator" analyzes trials generated
  outside the package (Weibull survival with a cure fraction) with
  `cutoff_fast()`, `analysis_fast()`, `switch_fast()`, and
  `simsummary_fast()`, and lists the columns each function reads.

* The group sequential design vignette is now evaluated when it is built
  instead of showing static output. The closed-form reference is computed with
  `gsDesign::gsProbability()` and the expected-event curve, and the
  non-proportional hazards example compares the log-rank and max-combo tests
  with `mc.alpha`. The rpact package is no longer used and is removed from
  Suggests.

* In the correlated PFS and OS vignette, the hierarchical rule now claims OS at
  the first look, from the PFS claim onward, at which OS crosses its boundary
  (Glimm, Maurer, and Bretz, 2010). The previous code used only the first OS
  crossing, so an OS crossing before the PFS claim prevented a later OS claim
  and the hierarchical OS power was understated. A duplicated row of the power
  table is removed.

* The multi-arm vignette computes the global-null error rates with
  `pairwise_fast()`, the multiregional vignette generates the regions from
  `dqrng` streams, and the Freidlin and Korn vignette draws the scenario with
  `gen_scenario_fast()` and uses `mc.alpha` for the max-combo test.

* The convergence statement of `coxph_fast()` (help page and README) is
  corrected. The residual error of the Pike-Halley Estimator is of the order
  of the cube of the error of the Pike anchor; it is negligible near the null
  hypothesis, but at a fixed hazard ratio away from 1 the Pike anchor keeps a
  bias that does not vanish with the sample size, so the residual error levels
  off at a small value rather than decreasing at the rate O_p(n^{-3/2})
  stated before.

* Help pages corrected or completed: the eligibility condition of
  `switch_fast()` for single-endpoint data includes dropout; the rules of
  `cutoff_fast()` for a look with only `max.time`, after an unreached look, and
  for looks out of order; the treatment of an unreached look in
  `analysis_fast()`, `pairwise_fast()`, and `simsummary_fast()`; the meaning of
  `intermediate` and `switched` in `simdata_fast()`; the scale of `std.err` in
  `survfit_fast()` compared with `survival::survfit()`; the `side` argument of
  `coxph_fast()`, whose return value has no p-value; the tests of
  `rmst_fast()` and `ahsw_fast()`, which follow `side`; and the description of
  `kmcurve_fast()` in the package overview, which has no risk table.

* The validation vignette compares `coxph_fast()` with `coxph()` using
  `ties = "breslow"`, the partial likelihood that the Pike-Halley Estimator
  approximates, and reports the difference from the computed values.

* The correlated PFS and OS vignette corrects the description of the events
  shared by the two endpoints and states that the OS information fractions are
  planning values and that the family-wise error rate is not evaluated there.

* Other vignette corrections: the Freidlin and Korn vignette describes the
  FH(0,1) weights correctly (bounded, but close to zero for early events); the
  multiregional vignette no longer attributes the hazard-reduction criterion
  to Teng et al. (2018), states the comparison of Method 1 and Method 2 as an
  observation, and removes a statement on scope that its three-region example
  contradicted; the multi-arm vignette reports the family-wise error rate with
  its Monte Carlo standard error and its normal approximation; the
  compare-logrank-rmst vignette draws the Kaplan-Meier curves of the data
  censored at the analysis; the speed comparison states the measurement
  conditions; and references listed but not cited are now cited in the text.

* README: the nominal levels of the workflow example are now those of the
  Lan-DeMets O'Brien-Fleming-type spending at 300 and 450 events; the return
  values of `milestone_fast()` and `ahr_fast()` are lists; and the agreement of
  `medsurv_fast(method = "nph")` with the nph package is stated for coinciding
  Kaplan-Meier and Nelson-Aalen medians.

## Tests

* New tests for `gen_scenario_fast()`, `plot.scenario_fast()`, and the print
  methods that had none, and for the `mc.alpha` shortcut of `analysis_fast()`.

* New tests for the fixes above: the variance of `wkm_fast()` for every weight
  against an independent reference and a calibration of its standard error by
  simulation; data with infinite latent times; switching and dropout in the
  illness-death model; degenerate milestone intervals; the order check of
  `presorted = TRUE`; and the input checks.

* New tests for subgroups in the illness-death model: a single cell
  reproduces the data without subgroups, each cell follows its own transition
  hazards, a subgroup matches its separate simulation with the accrual rate
  scaled by the prevalence, and the output works with `analysis_fast()`,
  `cutoff_fast()`, and `switch_fast()`.

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
