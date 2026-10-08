#' Fast Closed-Form Hazard Ratio Estimation via the Pike-Halley Estimator
#'
#' @description
#' Estimates the hazard ratio for a two-group parallel trial using the
#' Pike-Halley Estimator (Homma, 2025), a pure closed-form approximation to the Cox partial
#' likelihood maximizer. The function returns the point estimate, its standard
#' error on the log scale, and a Wald-type confidence interval, using output
#' names consistent with \code{summary(survival::coxph(...))}. The C++ backend
#' accepts pooled sorted vectors directly, performing group splitting and all
#' accumulation in a single C++ pass without intermediate R-level vector copies.
#'
#' @details
#' Let t_k (k = 1, ..., K) denote the distinct observed event times in the
#' pooled sample. At each t_k, let n_T_k and n_C_k be the numbers at risk
#' in the treatment and control groups just before t_k, and let O_T_k and
#' O_C_k be the numbers of events in each group, with n_k = n_T_k + n_C_k
#' and O_k = O_T_k + O_C_k. Define E_T = sum n_T_k O_k / n_k and
#' E_C = sum n_C_k O_k / n_k as the log-rank expected event totals.
#'
#' The Pike-Halley Estimator is obtained in three steps. First, the Pike
#' anchor is computed as theta_0 = (O_T E_C) / (O_C E_T). Second, the score
#' U_0, the observed information I_0, and the third-order curvature term J_0
#' of the Breslow partial likelihood are evaluated at eta_0 = log(theta_0):
#'
#' p_k = n_T_k theta_0 / (n_C_k + n_T_k theta_0)
#' U_0 = sum (O_T_k - O_k p_k)
#' I_0 = sum O_k p_k (1 - p_k)
#' J_0 = sum O_k p_k (1 - p_k) (1 - 2 p_k)
#'
#' Third, the closed-form Halley correction is applied:
#'
#' delta_hat = U_0 / I_0 - J_0 U_0^2 / (2 I_0^3)
#' theta_hat = theta_0 exp(delta_hat)
#'
#' The Halley step converges cubically, so the residual error
#' |log(theta_hat) - log(theta_Cox)| is of the order of the cube of the error
#' of the Pike anchor. Near the null hypothesis the anchor is
#' already close to the Cox estimate and the residual error is negligible. At
#' a fixed hazard ratio away from 1 the Pike anchor has a bias that does not
#' vanish as the sample size grows (Berry, Kitchin, and Mock, 1991), so the
#' residual error levels off at a small value that depends on the hazard ratio
#' and the censoring pattern instead of decreasing with the sample size.
#'
#' The Wald standard error on the log scale is SE = 1 / sqrt(I_0), where I_0
#' is the observed information evaluated at the Pike anchor. The Wald
#' confidence interval reported by \code{summary(coxph(...))} uses the
#' observed information at the maximum likelihood estimate instead. The two
#' differ by an amount of the order of the error of the Pike anchor, which is
#' small when the hazard ratio is not far from 1; Berry, Kitchin, and Mock
#' (1991) found the bias of the Pike estimator minimal for hazard ratios below
#' 3.
#'
#' With \code{strata}, the estimator targets the stratified Cox model, in which
#' each stratum has its own baseline hazard and the hazard ratio is common to
#' all strata. The at-risk sets are formed within each stratum, the observed
#' and expected totals O_T, O_C, E_T, and E_C are summed over strata to give
#' the stratified Pike anchor, and U_0, I_0, and J_0 are those of the
#' stratified Breslow partial likelihood (sums of the per-stratum terms). The
#' Halley correction and the Wald interval are then applied unchanged, so the
#' result approximates \code{coxph(Surv(time, event) ~ group + strata(s),
#' ties = "breslow")}. With a fixed number of strata the approximation error
#' shrinks at the same rate as in the unstratified case; it is larger when the
#' strata are small and the hazard ratio is far from 1, because the stratified
#' anchor is then further from the maximum likelihood estimate.
#'
#' The C++ core (\code{pihe_core}) accepts the pooled sorted data together
#' with an integer group indicator and performs group splitting, at-risk
#' counting, and per-distinct-event-time accumulation in a single left-to-right
#' scan. This eliminates the \code{rev(cumsum(rev(...)))}, \code{tapply()},
#' \code{which()}, \code{diff()}, and group-split vector copies present in the
#' pure-R version.
#'
#' The returned object has class \code{"coxph_fast"} and is a named numeric
#' vector of length 5. A \code{print()} method formats the result similarly
#' to \code{summary(coxph(...))}.
#'
#' @param time A numeric vector of follow-up times for all subjects (pooled
#'   over both groups).
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}. Any type
#'   that supports equality comparison is accepted.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group. Subjects with \code{group != control} are
#'   treated as the treatment group.
#' @param side 1 for a one-sided test in the direction of treatment benefit
#'   (hazard ratio below 1, i.e. a negative coefficient) or 2 for a two-sided
#'   test (default 2). The returned vector has no p-value: \code{side} is
#'   stored as an attribute, and the \code{print()} method reports the
#'   p-value that follows it. The confidence interval is always two-sided at
#'   \code{conf.level}.
#' @param conf.level A single numeric value in (0, 1) specifying the confidence
#'   level for the Wald interval. Defaults to 0.95.
#' @param presorted A logical value. If \code{TRUE}, \code{time},
#'   \code{event}, and \code{group} are assumed to be already sorted in
#'   ascending order of \code{time} (with \code{strata}, sorted by stratum
#'   and by time within stratum, so that the rows of each stratum are
#'   contiguous), and the internal \code{order()} call is skipped. If
#'   \code{FALSE} (default), sorting is handled internally.
#'   The order is checked, and an error is given when it does not hold.
#' @param strata An optional vector of stratum labels aligned with
#'   \code{time}, without missing values. When supplied, the stratified
#'   estimator described in Details is computed. Several stratification
#'   factors can be combined with \code{interaction()}. \code{NULL} (default)
#'   gives the unstratified estimator.
#'
#' @return An object of class \code{"coxph_fast"}, which is a named numeric
#'   vector of length 5 with elements matching the column names of
#'   \code{summary(coxph(...))$coefficients} and
#'   \code{summary(coxph(...))$conf.int}:
#' \describe{
#'   \item{\code{coef}}{Log hazard ratio log(theta_hat).}
#'   \item{\code{exp(coef)}}{Hazard ratio theta_hat (point estimate).}
#'   \item{\code{se(coef)}}{Standard error of \code{coef} on the log scale,
#'     equal to 1 / sqrt(I_0).}
#'   \item{\code{lower .95}}{Lower bound of the Wald confidence interval for
#'     the hazard ratio. The label reflects \code{conf.level} (e.g.,
#'     \code{"lower .90"} when \code{conf.level = 0.90}).}
#'   \item{\code{upper .95}}{Upper bound of the Wald confidence interval.}
#' }
#' With \code{strata}, the number of strata is stored in the attribute
#' \code{strata}.
#' Returns a vector of \code{NA_real_} values (still with class
#' \code{"coxph_fast"}) when the estimate cannot be computed (e.g., no
#' events, all events in one group, or \code{I_0 = 0}).
#'
#' @examples
#' library(survival)
#'
#' # Compare coxph_fast with coxph on the ovarian dataset.
#' # coxph() treats rx as numeric with rx=1 as the reference (control),
#' # so set control = 1 for a consistent comparison.
#' fit_fast <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
#' fit_fast
#'
#' fit_cox <- summary(coxph(Surv(futime, fustat) ~ rx, data = ovarian))
#' cat("coxph_fast HR :", fit_fast["exp(coef)"], "\n")
#' cat("coxph      HR :", fit_cox$coefficients[, "exp(coef)"], "\n")
#'
#' # Stratified by residual disease, compared with coxph(... + strata())
#' coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1,
#'            strata = ovarian$resid.ds)
#' coef(coxph(Surv(futime, fustat) ~ rx + strata(resid.ds), data = ovarian,
#'            ties = "breslow"))
#'
#' # presorted = TRUE: sort once outside, reuse inside a loop
#' ord <- order(ovarian$futime)
#' coxph_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
#'            control = 1, presorted = TRUE)
#'
#' \donttest{
#' # Speed comparison against coxph()
#' if (requireNamespace("microbenchmark", quietly = TRUE)) {
#'   microbenchmark::microbenchmark(
#'     coxph_fast = coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2),
#'     coxph      = coxph(Surv(futime, fustat) ~ rx, data = ovarian),
#'     times = 1000
#'   )
#' }
#' }
#'
#' @seealso
#' \code{\link[survival]{coxph}} for the standard iterative Cox estimator.
#' \code{\link{print.coxph_fast}} for the print method.
#'
#' @references
#' Cox, D. R. (1972). Regression models and life-tables. \emph{Journal of the
#' Royal Statistical Society. Series B (Methodological)}, \emph{34}(2),
#' 187-220.
#'
#' Berry, G., Kitchin, R. M., & Mock, P. A. (1991). A comparison of two simple
#' hazard ratio estimators based on the logrank test. \emph{Statistics in
#' Medicine}, \emph{10}(5), 749-755.
#'
#' Homma, G. (2025). One step from Pike to Cox: a closed-form hazard ratio
#' estimator. Manuscript under review.
#'
#' @importFrom stats qnorm setNames
#' @export
coxph_fast <- function(time, event, group, control, side = 2,
                       conf.level = 0.95, presorted = FALSE, strata = NULL) {

  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  if (length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be in (0, 1)")
  }

  # Prepare NA output with coxph-compatible names
  ci_lab <- conf.level * 100
  ci_lo  <- sprintf("lower .%g", ci_lab)
  ci_hi  <- sprintf("upper .%g", ci_lab)
  na_out <- setNames(
    rep(NA_real_, 5L),
    c("coef", "exp(coef)", "se(coef)", ci_lo, ci_hi)
  )

  # Input validation
  n <- length(time)
  if (length(event) != n || length(group) != n) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  use_strata <- !is.null(strata)
  if (use_strata) {
    if (length(strata) != n) {
      stop("'strata' must have the same length as 'time'")
    }
    if (anyNA(strata)) {
      stop("'strata' must not contain missing values")
    }
  }

  # Attach the attributes of a result (the number of strata when stratified)
  wrap <- function(v, n_strata = NULL) {
    out <- structure(v, conf.level = conf.level, side = side,
                     control = control, class = "coxph_fast")
    if (use_strata) attr(out, "strata") <- n_strata
    out
  }

  if (n == 0L) return(wrap(na_out, 0L))
  check_time_event(time, event)

  # Treatment indicator: 1 = treatment, 0 = control
  j <- two_group_indicator(group, control)

  # Stratum labels mapped to contiguous integers 1..S
  if (use_strata) {
    strata_int <- match(strata, sort(unique(strata)))
    n_strata   <- max(strata_int)
  } else {
    n_strata <- NULL
  }

  if (sum(event) == 0L) return(wrap(na_out, n_strata))

  # Sort by time (by stratum, then time, when stratified) when not presorted
  if (!presorted) {
    ord   <- if (use_strata) order(strata_int, time) else order(time)
    time  <- time[ord]
    event <- as.integer(event[ord])
    j     <- j[ord]
    if (use_strata) strata_int <- strata_int[ord]
  } else {
    check_presorted(time, if (use_strata) strata_int else NULL)
    event <- as.integer(event)
  }

  # C++ core: single scan over pooled sorted data (per stratum when
  # stratified) -> c(theta_0, U_0, I_0, J_0)
  res <- if (use_strata) {
    pihe_core_strat(time, event, j, strata_int)
  } else {
    pihe_core(time, event, j)
  }

  if (anyNA(res)) return(wrap(na_out, n_strata))

  theta_0 <- res[1L]
  U_0     <- res[2L]
  I_0     <- res[3L]
  J_0     <- res[4L]

  # Halley correction
  delta     <- U_0 / I_0 - (J_0 * U_0 * U_0) / (2 * I_0 * I_0 * I_0)
  theta_hat <- theta_0 * exp(delta)

  # Wald SE and CI on the log scale
  se_coef <- 1 / sqrt(I_0)
  coef    <- log(theta_hat)
  z       <- qnorm(1 - (1 - conf.level) / 2)

  out <- setNames(
    c(coef, theta_hat, se_coef,
      exp(coef - z * se_coef),
      exp(coef + z * se_coef)),
    c("coef", "exp(coef)", "se(coef)", ci_lo, ci_hi)
  )

  wrap(out, n_strata)
}
