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
