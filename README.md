# FastSurvival

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/FastSurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/FastSurvival/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

FastSurvival provides fast alternatives to the standard survival analysis
functions in the
[survival](https://cran.r-project.org/package=survival) package.
Every function is designed for repeated evaluation inside large simulation
loops — adaptive sample-size re-estimation, probability-of-success
calculations, regional consistency evaluation in multi-regional trials —
where the standard iterative or object-building overhead of `survfit()`,
`survdiff()`, and `coxph()` becomes a computational bottleneck. Core
computations are implemented in C++ via
[Rcpp](https://cran.r-project.org/package=Rcpp) for maximum performance.

## Functions

| Function | Replaces | Approximate speed gain |
|----------|----------|------------------------|
| `survfit_fast()` | `survfit()` + `summary()` at a single time point | ~50x |
| `survdiff_fast()` | `survdiff()` | ~40x |
| `coxph_fast()` | `coxph()` (point estimate + Wald CI) | ~30x |
| `simdata_fast()` | Custom simulation scripts | — |

Speed gains are based on median times from 1,000 microbenchmark replicates
on a typical phase-3 trial dataset (n = 600, event rate 80%) with
\code{presorted = TRUE}. Results vary by hardware and sample size.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("gosukehommaEX/FastSurvival")
```

## Quick start

```r
library(FastSurvival)
library(survival)

# ----------------------------------------------------------------
# survfit_fast(): Kaplan-Meier at a single time point
# ----------------------------------------------------------------
ord <- order(ovarian$futime)
t_s <- ovarian$futime[ord]
e_s <- ovarian$fustat[ord]

survfit_fast(t_s, e_s, t_eval = 500, conf.type = "log")

# ----------------------------------------------------------------
# survdiff_fast(): log-rank test
# ----------------------------------------------------------------
survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
              control = 1, side = 2)

# ----------------------------------------------------------------
# coxph_fast(): hazard ratio via the Pike-Halley Estimator (closed-form)
# ----------------------------------------------------------------
coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)

# ----------------------------------------------------------------
# simdata_fast(): clinical trial simulation
# ----------------------------------------------------------------

# Two-group trial, simple exponential, no dropout
df <- simdata_fast(
  nsim     = 1000,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.median = list(18, 24),
  seed     = 1
)
head(df)

# Two-group trial, piecewise exponential (delayed treatment effect)
df2 <- simdata_fast(
  nsim     = 1000,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.hazard = list(c(0.08, 0.08), c(0.08, 0.04)),
  e.time   = c(0, 6, Inf),
  seed     = 2
)
```

The three analysis functions return S3-class objects
(`survfit_fast`, `survdiff_fast`, `coxph_fast`) with `print()` methods that
display the results in a format similar to the corresponding `survival`
package functions. Each object is internally a numeric vector, so it can be
used directly in arithmetic, subsetting, and aggregation after stripping
the class (see the simulation example below).

## Design principles

**survfit_fast** evaluates the Kaplan-Meier estimator at a single specified
time point. A C++ backend locates the evaluation cutoff via binary search,
then accumulates the Kaplan-Meier product and the Greenwood variance sum in
a single scan over event positions, without constructing intermediate
vectors. This makes it orders of magnitude faster than `survfit()` plus
`summary()` when the same sorted data are evaluated repeatedly at a fixed
landmark time inside a simulation loop.

**survdiff_fast** computes the log-rank statistic using a two-pointer merge
scan over the pooled sorted data, walking the time axis once while
maintaining per-group at-risk counters and processing tied event times
atomically. The C++ backend avoids the rank construction, `tabulate()`, and
reverse cumulative sum operations of the standard R implementation, and the
function itself bypasses the S3 object construction and formula parsing of
`survdiff()`. It returns either a one-sided Z-score (`side = 1`) or a
two-sided chi-square statistic (`side = 2`).

**coxph_fast** implements the Pike-Halley Estimator proposed by Homma
(2025), a closed-form approximation to the Cox partial likelihood maximizer.
The estimator anchors at the Pike closed-form estimate and applies a single
analytic Halley correction to the Cox score, giving residual error of order
O_p(n^{-3/2}) relative to the Cox maximum likelihood estimate. On the
`pharmacoSmoking` dataset (tie rate 77.5%), the Pike-Halley Estimator
reproduces the Breslow-based Cox estimate to within on the order of 1e-08.
The Wald confidence interval uses the observed information at the Pike
anchor as the variance estimate. The C++ backend performs group splitting,
at-risk counting, and per-distinct-event-time accumulation in a single pass.

**simdata_fast** generates individual patient data for one- or two-group
time-to-event trials. Accrual times follow a piecewise uniform distribution.
Survival and dropout times follow either a simple or piecewise exponential
distribution, selected automatically based on whether a scalar or vector
hazard is supplied. C++ backends handle piecewise sampling and two-group
interleaving, and random number generation uses
[dqrng](https://cran.r-project.org/package=dqrng) for speed.

## Using the functions together

A typical simulation workflow combines all four functions. Inside the loop
the returned objects can be stored as-is. When aggregating results across
simulations, strip the S3 class with `unclass()` so that `rbind()` produces
an ordinary numeric matrix:

```r
library(FastSurvival)

df <- simdata_fast(
  nsim     = 1000,
  n        = c(100, 100),
  a.time   = c(0, 12),
  a.rate   = 200 / 12,
  e.hazard = list(0.08, 0.05),
  d.hazard = list(0.01, 0.01),
  seed     = 42
)

results <- vector("list", 1000L)
for (s in seq_len(1000L)) {
  d <- df[df$sim == s, ]
  results[[s]] <- coxph_fast(d$tte, d$event, d$group, control = 1)
}

hr_mat <- do.call(rbind, lapply(results, unclass))
colMeans(hr_mat)
```

## References

Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
estimator. *Manuscript under review.*

Collett, D. (2014). *Modelling Survival Data in Medical Research* (3rd ed.).
Chapman and Hall/CRC.

## License

MIT © 2025 Gosuke Homma
