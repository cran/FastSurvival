#' Compare Milestone Survival Probabilities Between Two Groups
#'
#' Compares the Kaplan-Meier survival probabilities of two groups at a
#' prespecified milestone timepoint. The point estimate of interest is the
#' difference in milestone survival, treatment minus control. Three inference
#' methods are provided. The \code{"wald"} method uses the unpooled Greenwood
#' variance directly. The \code{"loglog"} and \code{"mover"} methods build the
#' confidence interval for the difference with the method of variance estimates
#' recovery (MOVER), recovering the variance from the one-sample complementary
#' log-log and log transformed confidence intervals respectively. When the
#' estimate of a group is 0 or 1 at \code{tau}, its one-sample interval
#' degenerates to the estimate, so the interval of the difference is still
#' given; the \code{"loglog"} test statistic is then \code{NA}. See Tang
#' (2021) for the MOVER difference interval and Tang (2022) for the use of
#' milestone survival in trial design.
#'
#' The \code{"wald"} statistic is the difference divided by its Greenwood
#' standard error. The \code{"loglog"} statistic is the difference of the
#' complementary log-log transforms of the two estimates divided by its
#' standard error, so it is negative when the treatment group has the higher
#' milestone survival, and its one-sided p-value is the lower tail. The
#' \code{"mover"} statistic is the upper \code{(1 - conf.level) / 2} quantile
#' of the standard normal distribution (1.96 for \code{conf.level = 0.95})
#' times the difference divided by the half-width of the MOVER interval on the
#' side of zero. It is obtained by inverting the interval at
#' \code{conf.level}, so the statistic and its p-value depend on
#' \code{conf.level}, and the two-sided p-value is below
#' \code{1 - conf.level} exactly when the interval excludes zero. For
#' \code{"wald"} and \code{"mover"} a positive statistic favors treatment and
#' the one-sided p-value is the upper tail.
#'
#' @param time A numeric vector of follow-up times.
#'   Alternatively, a formula \code{Surv(time, event) ~ group}, with the
#'   variables in \code{data}; see \code{\link{formula_interface}}.
#' @param event An integer vector of event indicators, 1 for an event and 0
#'   for a censored observation.
#' @param group A vector with exactly two distinct values identifying the
#'   group.
#' @param control The value of \code{group} that denotes the control group. The
#'   other value is the treatment group and the difference is reported as
#'   treatment minus control.
#'   With the formula interface it defaults to the first level of the
#'   group variable.
#' @param side 1 for a one-sided test in the direction of treatment benefit
#'   (treatment milestone survival larger than control) or 2 for a two-sided
#'   test (default 2). The confidence interval is always reported as a
#'   two-sided interval at \code{conf.level}.
#' @param conf.level The confidence level for the reported intervals.
#' @param tau The milestone timepoint at which the survival probabilities are
#'   compared. A single positive number. A warning is given when it exceeds the
#'   largest observed time of either group, because the Kaplan-Meier curve is
#'   not estimated beyond that time.
#' @param method The inference method for the difference in milestone survival,
#'   one of \code{"wald"}, \code{"loglog"}, or \code{"mover"}.
#' @param presorted Logical. If \code{TRUE} the input is assumed to be sorted
#'   by \code{time} in ascending order and the internal sort is skipped. This
#'   is intended for repeated calls inside simulation loops.
#'   The order is checked, and an error is given when it does not hold.
#' @param data A data frame in which the variables of a formula given as
#'   \code{time} are evaluated. Used only with the formula interface; see
#'   \code{\link{formula_interface}}.
#'
#' @return An object of class \code{"milestone_fast"}, a list with the
#'   per-group milestone survival estimates and standard errors, the difference
#'   estimate with its confidence interval, the test statistic, and the p-value.
#'   The object also has the class \code{"fastsurv"}, whose
#'   \code{as.data.frame()}, \code{coef()}, and \code{confint()} methods
#'   are described in \code{\link{as.data.frame.fastsurv}}.
#'
#' @references
#' Tang, Y. (2021). Some new confidence intervals for Kaplan-Meier based
#' estimators from one and two sample survival data. \emph{Statistics in
#' Medicine}, \emph{40}(23), 4961-4976.
#'
#' Tang, Y. (2022). Complex survival trial design by the product integration
#' method. \emph{Statistics in Medicine}, \emph{41}(4), 798-814.
#'
#' @examples
#' set.seed(1)
#' time <- c(rexp(50, 0.1), rexp(50, 0.07))
#' event <- rep(1, 100)
#' group <- rep(c(0, 1), each = 50)
#' milestone_fast(time, event, group, control = 0, tau = 10, method = "loglog")
#'
#' # Formula interface
#' d <- data.frame(time, event, group)
#' milestone_fast(Surv(time, event) ~ group, data = d, tau = 10)
#'
#' @export
milestone_fast <- function(time, event, group, control, side = 2,
                           conf.level = 0.95, tau,
                           method = c("wald", "loglog", "mover"),
                           presorted = FALSE, data = NULL) {
  # Formula interface: Surv(time, event) ~ group evaluated in 'data'
  if (inherits(time, "formula")) {
    f <- fastsurv_formula(time, data,
                          event = if (!missing(event)) event,
                          group = if (!missing(group)) group,
                          group_mode = "required", fun = "milestone_fast()")
    if (missing(control) || is.null(control)) control <- f$control
    return(milestone_fast(f$time, f$event, f$group, control = control,
                          side = side, conf.level = conf.level, tau = tau,
                          method = method, presorted = presorted))
  }
  if (!is.null(data)) stop("'data' is used only with a formula")

  method <- match.arg(method)
  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided).")
  }

  if (!is.numeric(time)) stop("'time' must be numeric.")
  if (length(event) != length(time) || length(group) != length(time)) {
    stop("'time', 'event', and 'group' must have the same length.")
  }
  if (anyNA(time) || anyNA(group)) {
    stop("'time' and 'group' must not contain missing values.")
  }
  if (any(time < 0)) stop("'time' must be non-negative.")
  if (!(is.numeric(event) || is.logical(event))) {
    stop("'event' must be numeric or logical, not a factor or a character ",
         "vector.")
  }
  # Validate before converting, so that a value such as 0.7 is not truncated.
  if (anyNA(event) || !all(event == 0 | event == 1)) {
    stop("'event' must contain only 0 (censored) and 1 (event).")
  }
  event <- as.integer(event)
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("'tau' must be a single positive number.")
  }
  if (length(conf.level) != 1L || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number in (0, 1).")
  }

  # Map the two group levels to control (0) and treatment (1).
  if (is.factor(group)) {
    lev <- levels(droplevels(group))
  } else {
    lev <- sort(unique(group))
  }
  if (length(lev) != 2L) {
    stop("'group' must have exactly two distinct values.")
  }
  if (missing(control) || length(control) != 1L ||
      !(as.character(control) %in% as.character(lev))) {
    stop("'control' must be one of the two values in 'group'.")
  }
  trt <- lev[as.character(lev) != as.character(control)]
  grp01 <- as.integer(as.character(group) != as.character(control))
  check_tau_follow_up(as.numeric(time), grp01, tau)
  if (presorted) check_presorted(time)

  core <- milestone_core(as.numeric(time), event, grp01,
                         as.numeric(tau), as.logical(presorted))

  surv0 <- core$surv0
  surv1 <- core$surv1
  v0 <- core$var0
  v1 <- core$var1
  se0 <- sqrt(v0)
  se1 <- sqrt(v1)

  z <- stats::qnorm(1 - (1 - conf.level) / 2)

  # One-sample confidence limits for each group under the chosen transform.
  one_sample_ci <- function(surv, se) {
    if (is.na(surv) || is.na(se)) return(c(NA_real_, NA_real_))
    if (method == "wald") {
      return(c(surv - z * se, surv + z * se))
    }
    # An estimate of 0 or 1 has zero Greenwood variance, and its interval
    # degenerates to the estimate.
    if (surv <= 0 || surv >= 1) return(c(surv, surv))
    if (method == "mover") {
      # log transform: var(log surv) = se^2 / surv^2
      g_se <- se / surv
      return(c(surv * exp(-z * g_se), surv * exp(z * g_se)))
    }
    # method == "loglog": g = log(-log(surv))
    g_se <- se / (surv * abs(log(surv)))
    return(c(surv^exp(z * g_se), surv^exp(-z * g_se)))
  }

  ci0 <- one_sample_ci(surv0, se0)
  ci1 <- one_sample_ci(surv1, se1)

  diff_est <- surv1 - surv0

  if (method == "wald") {
    se_diff <- sqrt(v0 + v1)
    diff_lower <- diff_est - z * se_diff
    diff_upper <- diff_est + z * se_diff
    stat <- diff_est / se_diff
  } else {
    # MOVER difference interval (Tang 2021, equation 10).
    l0 <- ci0[1]; u0 <- ci0[2]; l1 <- ci1[1]; u1 <- ci1[2]
    sigma_l <- (l1 - surv1)^2 + (u0 - surv0)^2
    sigma_u <- (u1 - surv1)^2 + (l0 - surv0)^2
    diff_lower <- diff_est - sqrt(sigma_l)
    diff_upper <- diff_est + sqrt(sigma_u)
    if (method == "loglog") {
      # Variance-stabilised complementary log-log test of equal survival.
      if (surv0 > 0 && surv0 < 1 && surv1 > 0 && surv1 < 1) {
        g0 <- log(-log(surv0))
        g1 <- log(-log(surv1))
        vg0 <- v0 / (surv0 * log(surv0))^2
        vg1 <- v1 / (surv1 * log(surv1))^2
        stat <- (g1 - g0) / sqrt(vg0 + vg1)
      } else {
        stat <- NA_real_
      }
    } else {
      # MOVER binding-side test statistic: the interval at conf.level inverted,
      # so the statistic depends on conf.level.
      bind <- if (diff_est >= 0) sigma_l else sigma_u
      stat <- if (is.finite(bind) && bind > 0) z * diff_est / sqrt(bind) else NA_real_
    }
  }

  # One-sided p-value (side = 1) tests the treatment-benefit direction, i.e.
  # higher milestone survival in the treatment group. The benefit tail depends
  # on the sign convention of the test statistic: for "wald" and "mover" a
  # positive statistic favours treatment (upper tail), whereas for "loglog" the
  # complementary log-log transform reverses the sign so benefit lies in the
  # lower tail.
  p_value <- if (is.na(stat)) {
    NA_real_
  } else if (side == 2L) {
    2 * stats::pnorm(-abs(stat))
  } else if (method == "loglog") {
    stats::pnorm(stat)
  } else {
    stats::pnorm(-stat)
  }

  structure(
    list(
      surv = c(control = surv0, treatment = surv1),
      std.err = c(control = se0, treatment = se1),
      surv.lower = c(control = ci0[1], treatment = ci1[1]),
      surv.upper = c(control = ci0[2], treatment = ci1[2]),
      n = c(control = core$n0, treatment = core$n1),
      n.event = c(control = core$n_event0, treatment = core$n_event1),
      n.risk = c(control = core$n_risk0, treatment = core$n_risk1),
      diff = diff_est,
      diff.lower = diff_lower,
      diff.upper = diff_upper,
      statistic = stat,
      p.value = p_value,
      tau = as.numeric(tau),
      method = method,
      side = side,
      conf.level = conf.level,
      group.labels = c(control = as.character(control),
                       treatment = as.character(trt))
    ),
    class = c("milestone_fast", "fastsurv")
  )
}
