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
