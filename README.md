# FastSurvival <img src="man/figures/logo.png" align="right" height="139" alt="FastSurvival logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/FastSurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/FastSurvival/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/FastSurvival)](https://CRAN.R-project.org/package=FastSurvival)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/FastSurvival)](https://cran.r-project.org/package=FastSurvival)
[![CRAN downloads total](https://cranlogs.r-pkg.org/badges/grand-total/FastSurvival)](https://cran.r-project.org/package=FastSurvival)
[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://gosukehommaEX.github.io/FastSurvival/)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/FastSurvival/graph/badge.svg)](https://app.codecov.io/gh/gosukehommaEX/FastSurvival)
<!-- badges: end -->

FastSurvival provides fast alternatives to the standard survival analysis
functions in the
[survival](https://cran.r-project.org/package=survival) package, together
with a simulation layer for designing and evaluating time-to-event trials.
Every function is designed for repeated evaluation inside large simulation
loops — adaptive sample-size re-estimation, probability-of-success
calculations, regional consistency evaluation in multi-regional trials —
where the standard iterative or object-building overhead of `survfit()`,
`survdiff()`, and `coxph()` becomes a computational bottleneck. Core
computations are implemented in C++ via
[Rcpp](https://cran.r-project.org/package=Rcpp) for maximum performance.
Collett (2014) gives the methodological background of the standard methods.

The full documentation site, including function references and articles, is
available at <https://gosukehommaEX.github.io/FastSurvival/>.

## Functions

### Estimation and testing

| Function | Description |
|----------|-------------|
| `survfit_fast()` | Single-time-point Kaplan-Meier estimate with Greenwood standard error and plain / log / log-log confidence intervals. |
| `survdiff_fast()` | Log-rank test, including weighted (Fleming-Harrington, modestly-weighted, Gehan-Breslow, Tarone-Ware) and stratified variants. |
| `coxph_fast()` | Closed-form hazard ratio via the Pike-Halley Estimator with a Wald confidence interval, optionally stratified. |
| `rmst_fast()` | Restricted mean survival time for a single group or a two-group comparison (difference and ratio). |
| `wmst_fast()` | Window mean survival time over an interval, generalizing the restricted mean survival time, for a single group or a two-group comparison. |
| `milestone_fast()` | Two-group comparison of Kaplan-Meier survival at a milestone timepoint (Wald, log-log, and MOVER methods). |
| `medsurv_fast()` | Median survival time for a single group or a two-group comparison, with a native and an `nph`-compatible variance method. |
| `maxcombo_fast()` | Max-combo test, the maximum over a set of Fleming-Harrington weighted log-rank statistics, for non-proportional hazards. |
| `rmw_fast()` | Robust modestly-weighted test (Magirr and Öhrn): the maximum of the standard log-rank and a modestly-weighted log-rank statistic, for non-proportional hazards. |
| `wkm_fast()` | Weighted Kaplan-Meier (Pepe-Fleming) test, the weighted integrated difference between two Kaplan-Meier curves, for non-proportional hazards. |
| `ahsw_fast()` | Average hazard with survival weight of Uno and Horiguchi, with ratio and difference contrasts. |
| `ahr_fast()` | Kalbfleisch-Prentice average hazard ratio between two groups over a restricted interval, for non-proportional hazards. |

### Simulation

| Function | Description |
|----------|-------------|
| `simdata_fast()` | Individual patient data simulator for one-, two-, or multi-arm trials, with piecewise-uniform accrual, piecewise-exponential survival and dropout, optional subgroups, and two correlated endpoints from an illness-death model, which can also be combined with subgroups. |
| `cutoff_fast()` | Per-simulation calendar cutoffs of the analysis looks from combined trigger rules (target events, planned and maximum calendar times, minimum time after the previous look, minimum follow-up after a number of enrolled subjects). |
| `analysis_fast()` | Interim or sequential analysis of simulated data at one or more looks, defined by target event counts, calendar times, or per-simulation cutoffs from `cutoff_fast()`. |
| `switch_fast()` | Treatment switching in simulated data: switching at an intermediate event (such as progression), at an opening time after an interim analysis, or at the later of the two, optionally only in the trials selected by an interim decision, with an accelerated-failure-time or new-hazard effect after the switch. |
| `pairwise_fast()` | Pairwise comparison of each experimental arm against a shared control on multi-arm data, at fixed calendar looks, a primary contrast's event-driven cutoffs, or supplied per-simulation cutoffs, with an optional Bonferroni adjustment. |
| `simsummary_fast()` | Operating-characteristic summary (rejection and futility rates, stopping-look distribution, expected timing) from `analysis_fast()` output and supplied boundaries. |

### Visualization

| Function | Description |
|----------|-------------|
| `gen_scenario_fast()` | Assemble one or more two-group scenarios into a `scenario_fast` object for design-stage exploration, with `plot()` (analytic survival curves and the piecewise hazard ratio) and `print()` methods. |
| `kmcurve_fast()` | Build the Kaplan-Meier curves of a single trial realization into a `kmcurve_fast` object, with `plot()` (curves, restricted-mean shading, and a smoothed time-varying hazard-ratio panel) and `print()` methods. |

The estimation and testing functions return S3-class objects with `print()`
methods that display the results in a format similar to the corresponding
`survival` package functions.

## Installation

```r
# Install the released version from CRAN
install.packages("FastSurvival")

# Or install the development version from GitHub
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
# rmst_fast(): restricted mean survival time (two-group)
# ----------------------------------------------------------------
rmst_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
          control = 1, tau = 500)

# ----------------------------------------------------------------
# wmst_fast(): window mean survival time over [100, 500]
# ----------------------------------------------------------------
wmst_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
          control = 1, tau1 = 100, tau2 = 500)

# ----------------------------------------------------------------
# medsurv_fast(): median survival time (two-group)
# ----------------------------------------------------------------
medsurv_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)

# ----------------------------------------------------------------
# wkm_fast(): weighted Kaplan-Meier (Pepe-Fleming) test
# ----------------------------------------------------------------
wkm_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)

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

# Two correlated endpoints (e.g. PFS and OS) from an illness-death model:
# an intermediate event (h01) and a terminal event (h02). With the post-event
# hazard defaulting to h02 this is the Fleischer maximal-independence model.
df3 <- simdata_fast(
  nsim       = 1000,
  n          = c(100, 100),
  a.time     = c(0, 12),
  a.rate     = 200 / 12,
  h01.median = list(6, 9),     # intermediate event (e.g. progression)
  h02.median = list(14, 18),   # terminal event (e.g. death)
  seed       = 3
)
```

## Design principles

**survfit_fast** evaluates the Kaplan-Meier estimator at a single specified
time point. A C++ backend locates the evaluation cutoff via binary search,
then accumulates the Kaplan-Meier product and the Greenwood variance sum in
a single scan over event positions, without constructing intermediate
vectors. This makes it orders of magnitude faster than `survfit()` plus
`summary()` when the same sorted data are evaluated repeatedly at a fixed
landmark time inside a simulation loop.

**survdiff_fast** computes the log-rank statistic (Mantel, 1966) using a
two-pointer merge
scan over the pooled sorted data, walking the time axis once while
maintaining per-group at-risk counters and processing tied event times
atomically. The C++ backend avoids the rank construction, `tabulate()`, and
reverse cumulative sum operations of the standard R implementation, and the
function itself bypasses the S3 object construction and formula parsing of
`survdiff()`. Beyond the ordinary unweighted test, it supports the
Fleming-Harrington, modestly-weighted, Gehan-Breslow, and Tarone-Ware
weighted tests, as well as stratified and stratified-weighted variants, all
sharing the same single-scan backend.

**coxph_fast** implements the Pike-Halley Estimator proposed by Homma
(2025), a closed-form approximation to the Cox partial likelihood maximizer.
The estimator anchors at the Pike closed-form estimate, the ratio of the
observed-to-expected event ratios of the two groups from the log-rank
computation (Berry, Kitchin, and Mock, 1991), and applies a single analytic
Halley correction to the Cox score. Because the Halley step converges
cubically, the residual error relative to the Cox maximum likelihood estimate
is of the order of the cube of the error of the Pike anchor. It is negligible
near the null hypothesis. At a fixed hazard ratio away from 1 the Pike anchor
keeps a bias that does not vanish with the sample size, so the residual error
levels off at a small value. On the
`pharmacoSmoking` dataset (tie rate 77.5%), the Pike-Halley Estimator
reproduces the Breslow-based Cox estimate to within on the order of 1e-08.
The Wald confidence interval uses the observed information at the Pike
anchor as the variance estimate. The C++ backend performs group splitting,
at-risk counting, and per-distinct-event-time accumulation in a single pass.
With `strata`, the risk sets are formed within each stratum and the Pike
anchor, score, and information are summed over strata, which approximates the
stratified Cox model `coxph(... + strata(s), ties = "breslow")`.

**rmst_fast** integrates the Kaplan-Meier survival step function up to a
horizon (Royston and Parmar, 2013) in a single C++ scan, reused once per
group. With a single group it
returns the RMST with a Greenwood-type standard error and a Wald confidence
interval; with two groups it adds the RMST difference and ratio contrasts,
each with a standard error, confidence interval, and test.

**wmst_fast** computes the window mean survival time, the area under the
Kaplan-Meier curve between a lower and an upper window limit, generalizing the
restricted mean survival time (which is the special case with a lower limit of
zero). It targets late differences and early crossings by focusing the mean
on a clinically relevant window, following Paukner and Chappell (2021). The
same single C++ scan as `rmst_fast` accumulates the windowed area and a
Greenwood-type variance in which each event time contributes its squared
remaining window area. With two groups it reports the difference (treatment
minus control) with a standard error, confidence interval, and test.

**milestone_fast** compares the Kaplan-Meier survival probabilities of two
groups at a prespecified milestone timepoint, estimating the difference
(treatment minus control). Three inference methods are provided: a Greenwood
Wald interval, and two methods that build the difference interval by the
method of variance estimates recovery (MOVER) from complementary log-log and
log transformed one-sample intervals, following Tang (2021).

**medsurv_fast** estimates the Kaplan-Meier median survival time and, for two
groups, the difference (treatment minus control), in a single C++ scan that
locates the median on the Kaplan-Meier curve. Because the variance of a
quantile depends on the local hazard, two variance methods are provided: a
native method using a kernel-smoothed hazard with a Greenwood increment, and an
`nph`-compatible method using a local constant hazard, which reproduces the
standard error and p-value of the median comparison in the `nph` package to
numerical precision when the Kaplan-Meier and Nelson-Aalen medians coincide
(the usual case). The point estimate is the same under both methods.

**maxcombo_fast** computes the max-combo test, the maximum over a set of
Fleming-Harrington weighted log-rank statistics, for comparing survival
under non-proportional hazards (Karrison, 2016; Lin et al., 2020). The C++
backend evaluates every weighted
numerator and the full between-scheme covariance matrix in a single scan,
and the p-value follows from the implied multivariate normal distribution.
Taking the most extreme of several complementary weights makes the test
robust to the shape of the hazard difference, with the multiplicity
accounted for through the joint distribution.

**rmw_fast** computes the robust modestly-weighted (rMW) test of Magirr and
Öhrn (2026), the maximum of the standard log-rank statistic and a single
modestly-weighted log-rank statistic. The modestly-weighted component uses
the survival-threshold parameterization, with the weight capped at the
reciprocal of `s_star` (so `s_star = 0.5` caps the weight at 2 and `s_star = 1`
recovers the log-rank test). The C++ backend computes both component
numerators, both variances, and their null covariance in a single scan, and
the combined statistic and p-value follow from the implied bivariate normal
distribution. Because the standard log-rank statistic is one of the two
components and the components are strongly correlated under the null, the
multiplicity adjustment is small, so the test loses little power relative to
the log-rank test in worst-case scenarios while gaining power under delayed
effects.

**wkm_fast** computes the weighted Kaplan-Meier test of Pepe and Fleming
(1989), the weighted integral of the difference between the two Kaplan-Meier
curves, standardized to a Wald z statistic. Unlike weighted log-rank tests it
targets the integrated difference in survival and remains sensitive when the
hazards cross. The default Pepe-Fleming combined weight, built from the
censoring distributions of the two groups, gives less weight to regions with
heavy
censoring; square-root and constant weights are also available. The single C++
scan evaluates the per-group Kaplan-Meier curves, the censoring weights, and
the variance over the pooled event grid, and the result reproduces the
weighted Kaplan-Meier statistic of the `nphsim` package.

**ahsw_fast** computes the average hazard with survival weight of Uno and
Horiguchi for two groups. The average hazard on the window up to a horizon
is the ratio of the cumulative event probability to the restricted mean
survival time, both from the Kaplan-Meier estimate. The function returns the
per-group average hazard, the ratio of average hazards (RAH) on the log
scale, and the difference of average hazards (DAH) on the identity scale,
each with a confidence interval and a two-sided test.

**ahr_fast** computes the average hazard ratio of Kalbfleisch and Prentice
(1981) between two groups over a restricted interval, the estimator used by
Dormuth et al. (2024) for sample-size calculation under non-proportional
hazards. Writing the two survival functions as S0 and S1, each group's share
of the total hazard is obtained by integrating one survival curve against the
increments of the other; the average hazard ratio is the ratio of the two
shares, and equals the true hazard ratio under proportional hazards. A C++
backend builds both Kaplan-Meier curves and their Greenwood variances on the
pooled event grid, then evaluates the shares, their direct variances, and the
double-integral covariance terms in linear scans. The primary test is on the
group-share scale, as in the source method, with an equivalent test and a
confidence interval reported on the log scale.

**simdata_fast** generates individual patient data for one-group, two-group,
or multi-arm time-to-event trials. Accrual times follow a piecewise-uniform distribution,
and survival and dropout times follow either a simple or piecewise
exponential distribution, selected automatically based on whether a scalar
or vector hazard is supplied. Optional subgroups are defined by a prevalence
specification. It can also generate two correlated time-to-event
endpoints, such as progression-free and overall survival, from an
illness-death model with transition hazards for the intermediate event,
the direct terminal event, and the post-event terminal event, and with
optional treatment switching. The entire generation pipeline runs in a
single C++ kernel
that materializes the output data frame once, avoiding intermediate R-level
vector operations and copies, and random-number generation uses
[dqrng](https://cran.r-project.org/package=dqrng) for speed. With the
`stream` argument a large study can be generated in batches from independent
random-number streams, sequentially or in parallel, with results that do not
depend on how the batches are distributed.

**cutoff_fast** determines the calendar time of every analysis look in every
simulated trial in a single C++ pass, from a target number of events (on the
whole trial, a subset of subjects, or another endpoint), a planned calendar
time, a maximum calendar time, a minimum time after the previous look, and a
minimum follow-up after a given number of enrolled subjects, combined as in
`get_analysis_date()` of the
[simtrial](https://cran.r-project.org/package=simtrial) package. The resulting
cutoffs are passed to `analysis_fast()`, `pairwise_fast()`, and
`switch_fast()`, so that the same analysis times can drive several endpoints,
several contrasts, and a modification of the data after an interim analysis.

**analysis_fast** performs interim or sequential analyses of simulated
two-group data at one or more looks, defined by a target cumulative number
of events, by a calendar time, or by per-simulation cutoffs from
`cutoff_fast()`. At every look the data are
administratively censored at the corresponding calendar cutoff, and the
requested statistics are computed for each simulated trial. The censoring,
time sorting, and per-cell statistics are handled by a fused C++ kernel that
reuses the same analysis cores as the standalone functions, so the same
results are obtained without the per-iteration overhead of repeated wrapper
calls. Statistics can also be reported within each subgroup. The max-combo
p-value, a multivariate normal integral per trial and look, dominates the
computing time; when only the decisions at nominal levels are needed,
`mc.alpha` integrates only the p-values that the Bonferroni bounds do not
already place on one side of the level, with unchanged decisions.

**pairwise_fast** compares each experimental arm against a shared control on
multi-arm data by running `analysis_fast()` once per contrast and stacking the
results with an `arm` column. With fixed calendar looks every contrast shares
the same data cutoff. With an event-driven design a designated primary contrast
fixes the per-simulation calendar cutoff at which all contrasts are analyzed,
reproducing the standard rule that the primary analysis defines a single data
cutoff. Per-simulation cutoffs from `cutoff_fast()` can also be supplied
directly. An optional Bonferroni adjustment controls the family-wise error rate
across the contrasts at each look.

**switch_fast** applies a treatment-switching rule to simulated data by
changing only the outcomes after each subject's switch: switching at the
intermediate event of an illness-death simulation (for example at
progression), at an opening time after an interim analysis (crossover at a
milestone), or at the later of the two, with the remaining time to the terminal
event multiplied by an acceleration factor (the causal model of the
rank-preserving structural failure time method of Robins and Tsiatis, 1991) or
redrawn from a new hazard. Because
the history before the switch is kept, an interim analysis is unaffected, so a
crossover decided by an interim result is simulated for all trials at once:
analyze the interim, select the trials, switch, and analyze the later looks.

**simsummary_fast** aggregates the per-simulation, per-look output of
`analysis_fast()` into operating characteristics: the rejection rate, the
futility-stopping rate, the distribution of the stopping look, and the
expected analysis timing. It consumes group-sequential boundaries computed
elsewhere (for example by gsDesign or rpact) and estimates the stopping
probabilities by Monte Carlo over the simulated trials; it does not compute
or spend alpha itself.

**gen_scenario_fast** assembles one or more two-group scenarios, each a
complete set of `simdata_fast()` arguments, into a single object for
design-stage exploration. Its `plot()` method draws, for every scenario, the
analytic control and experimental survival curves with the piecewise-constant
hazard ratio on a secondary axis, computed directly from the
piecewise-exponential parameters without simulation, and its `print()` method
summarizes the medians, the hazard ratio at the start and end of the window,
and whether the curves cross. The same object can be passed to
`simdata_fast()` to generate the data.

**kmcurve_fast** builds the Kaplan-Meier step functions of a single realized
data set, for example one replicate from `simdata_fast()`, for both the
two-group and single-group cases. Its `plot()` method draws the curves and can
shade the area under each curve up to a truncation time to visualize the
restricted mean survival time, and can add a second panel showing the smoothed
time-varying hazard ratio, a kernel-smoothed Nelson-Aalen estimate, together
with the constant Cox proportional-hazards estimate. Its `print()` method
reports the sample sizes, events, and medians, with the Cox hazard ratio in
the two-group case.

## Using the functions together

A typical simulation workflow combines the simulation and analysis layers.
`simdata_fast()` generates the data once, `analysis_fast()` performs the
looks and computes the statistics for every simulated trial, and
`simsummary_fast()` aggregates the operating characteristics against the
supplied boundaries:

```r
library(FastSurvival)

# 1. Simulate a two-group trial with a delayed treatment effect
df <- simdata_fast(
  nsim     = 1000,
  n        = c(350, 350),
  a.time   = c(0, 12),
  a.rate   = 700 / 12,
  e.hazard = list(c(0.08, 0.08), c(0.08, 0.04)),
  e.time   = c(0, 6, Inf),
  d.hazard = list(0.01, 0.01),
  seed     = 42
)

# 2. Analyze at two interim looks defined by target event counts
res <- analysis_fast(
  df, control = 1,
  event.looks = c(300, 450),
  stat = "logrank", side = 1
)

# 3. Summarize operating characteristics against efficacy boundaries:
#    nominal one-sided levels of Lan-DeMets O'Brien-Fleming-type spending
#    (one-sided alpha 0.025) at 300 and 450 events, from
#    gsDesign::gsDesign(k = 2, test.type = 1, timing = 300 / 450,
#                       sfu = gsDesign::sfLDOF)
oc <- simsummary_fast(
  res,
  p.col = "logrank.p",
  alpha = c(0.0060, 0.0231)
)
oc
```

The estimation functions can also be used directly inside a loop. Except for
`milestone_fast()` and `ahr_fast()`, which return lists, each returned object
is internally a numeric vector, so after stripping the class with
`unclass()` the per-simulation results can be combined with `rbind()` into an
ordinary numeric matrix:

```r
results <- vector("list", 1000L)
for (s in seq_len(1000L)) {
  d <- df[df$sim == s, ]
  results[[s]] <- coxph_fast(d$tte, d$event, d$group, control = 1)
}

hr_mat <- do.call(rbind, lapply(results, unclass))
colMeans(hr_mat)
```

## References
 
Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete
observations. *Journal of the American Statistical Association*, 53(282),
457-481.
 
Mantel, N. (1966). Evaluation of survival data and two new rank order
statistics arising in its consideration. *Cancer Chemotherapy Reports*,
50(3), 163-170.
 
Gehan, E. A. (1965). A generalized Wilcoxon test for comparing arbitrarily
singly-censored samples. *Biometrika*, 52, 203-223.
 
Tarone, R. E., & Ware, J. (1977). On distribution-free tests for equality of
survival distributions. *Biometrika*, 64, 156-160.
 
Fleming, T. R., & Harrington, D. P. (1991). *Counting Processes and Survival
Analysis*. New York: John Wiley & Sons.
 
Pepe, M. S., & Fleming, T. R. (1989). Weighted Kaplan-Meier statistics: a
class of distance tests for censored survival data. *Biometrics*, 45(2),
497-507.
 
Pepe, M. S., & Fleming, T. R. (1991). Weighted Kaplan-Meier statistics: large
sample and optimality considerations. *Journal of the Royal Statistical
Society. Series B (Methodological)*, 53(2), 341-352.
 
Magirr, D., & Burman, C.-F. (2019). Modestly weighted logrank tests.
*Statistics in Medicine*, 38(20), 3782-3790.
 
Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal
Statistical Society. Series B (Methodological)*, 34(2), 187-220.
 
Berry, G., Kitchin, R. M., & Mock, P. A. (1991). A comparison of two simple
hazard ratio estimators based on the logrank test. *Statistics in Medicine*,
10(5), 749-755.
 
Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
estimator. *Manuscript under review.*
 
Royston, P., & Parmar, M. K. B. (2013). Restricted mean survival time: an
alternative to the hazard ratio for the design and analysis of randomized
trials with a time-to-event outcome. *BMC Medical Research Methodology*, 13,
152.
 
Paukner, M., & Chappell, R. (2021). Window mean survival time. *Statistics in
Medicine*, 40(25), 5521-5533.
 
Uno, H., Claggett, B., Tian, L., et al. (2014). Moving beyond the hazard
ratio in quantifying the between-group difference in survival analysis.
*Journal of Clinical Oncology*, 32(22), 2380-2385.
 
Tang, Y. (2021). Some new confidence intervals for Kaplan-Meier based
estimators from one and two sample survival data. *Statistics in Medicine*,
40(23), 4961-4976.
 
Tang, Y. (2022). Complex survival trial design by the product integration
method. *Statistics in Medicine*, 41(4), 798-814.
 
Karrison, T. G. (2016). Versatile tests for comparing survival curves based
on weighted log-rank statistics. *The Stata Journal*, 16(3), 678-690.
 
Lin, R. S., Lin, J., Roychoudhury, S., et al. (2020). Alternative analysis
methods for time to event endpoints under nonproportional hazards: a
comparative analysis. *Statistics in Biopharmaceutical Research*, 12(2),
187-198.
 
Magirr, D., & Öhrn, F. (2026). Robust modestly weighted log-rank tests.
*Pharmaceutical Statistics*, 25(1), e70066.
 
Uno, H., & Horiguchi, M. (2023). Ratio and difference of average hazard with
survival weight: new measures to quantify survival benefit of new therapy.
*Statistics in Medicine*, 42(7), 936-952.
 
Kalbfleisch, J. D., & Prentice, R. L. (1981). Estimation of the average
hazard ratio. *Biometrika*, 68(1), 105-112.
 
Dormuth, I., Pauly, M., Rauch, G., & Herrmann, C. (2024). Sample size
calculation under nonproportional hazards using average hazard ratios.
*Biometrical Journal*, 66(6), e202300271.
 
Robins, J. M., & Tsiatis, A. A. (1991). Correcting for non-compliance in
randomized trials using rank preserving structural failure time models.
*Communications in Statistics - Theory and Methods*, 20(8), 2609-2631.
 
Collett, D. (2014). *Modelling Survival Data in Medical Research* (3rd ed.).
Chapman and Hall/CRC.

## License

MIT © 2026 Gosuke Homma
