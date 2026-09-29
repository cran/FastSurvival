#' Fast non-parametric median survival time and between-group difference
#'
#' Estimates the Kaplan-Meier median survival time for a single group, or the
#' difference in median survival time between a treatment group and a control
#' group, together with standard errors, confidence intervals and a Wald test
#' for the difference. The estimator is designed for repeated evaluation inside
#' simulation loops, with the single scan over the sorted data performed in
#' C++.
#'
#' The median in each group is the first event time at which the Kaplan-Meier
#' estimate drops to 0.5 or below. When the estimate equals 0.5 exactly over an
#' interval, the midpoint of that interval is returned (the next event time, or
#' the last observed time when no later event exists, closes the interval). This
#' is the convention used by \code{survfit}, including its tolerance of
#' \code{sqrt(.Machine$double.eps)} when comparing the estimate with 0.5. The
#' point estimate is the same for both variance methods.
#'
#' Two variance methods are available through the \code{method} argument. With
#' \code{method = "km"} the variance of the estimated median follows the
#' counting process delta method, var(median) = greenwood_sum / hazard^2, where
#' the instantaneous hazard at the median is obtained by a Ramlau-Hansen kernel
#' smoother with an Epanechnikov kernel and bandwidth \code{bw}. With
#' \code{method = "nph"} the variance reproduces the computation used by
#' \code{nph::nphparams} with \code{param_type = "Q"} and
#' \code{haz_method = "local"}: a local constant hazard at the median and a
#' \code{sum 1 / (Y - k)^2} variance increment. The two methods give the same
#' median but generally different standard errors, since the variance of a
#' quantile depends on the local hazard estimate. When the Kaplan-Meier and
#' Nelson-Aalen medians coincide, which is the usual case, \code{method = "nph"}
#' reproduces the \code{nph::nphparams} standard error and p-value to numerical
#' precision.
#'
#' The difference is computed as treatment minus control, so a positive
#' difference indicates a longer median survival time under treatment.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Integer vector, 1 for an event and 0 for censoring.
#' @param group Optional grouping vector with exactly two distinct levels for a
#'   two-group comparison. If omitted, a single-group median is returned.
#' @param control The level of \code{group} that denotes the control group.
#'   Required when \code{group} is supplied.
#' @param side Either 2 for a two-sided test or 1 for a one-sided test of
#'   treatment superiority (difference greater than 0).
#' @param conf.level Confidence level for the intervals.
#' @param conf.type Confidence interval type for each group median, either
#'   \code{"plain"} or \code{"log"}. The interval for the difference is always
#'   on the plain scale.
#' @param method Variance method, either \code{"km"} (Greenwood increment with a
#'   kernel hazard) or \code{"nph"} (local constant hazard with a
#'   \code{sum 1 / (Y - k)^2} increment, matching \code{nph::nphparams}).
#' @param bw Optional kernel bandwidth for the hazard at the median, used only
#'   when \code{method = "km"}. Either a single value applied to both groups or
#'   one value per group. If omitted, a Silverman type default is used.
#' @param presorted Logical; set to \code{TRUE} when \code{time}, \code{event}
#'   and \code{group} are already sorted in ascending order of \code{time}, to
#'   skip the internal sort. Defaults to \code{FALSE}.
#'
#' @return A named numeric vector of class \code{"medsurv_fast"}. For a single
#'   group the elements are the median, its standard error and confidence
#'   limits. For two groups the elements are the control and treatment medians,
#'   their difference, the standard errors, the confidence limits and the Wald
#'   statistics for the difference.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' g <- rep(0:1, each = n / 2)
#' tt <- rexp(n, rate = ifelse(g == 0, 0.1, 0.07))
#' cc <- rexp(n, rate = 0.02)
#' time <- pmin(tt, cc)
#' event <- as.integer(tt <= cc)
#' medsurv_fast(time, event, group = g, control = 0)
#' medsurv_fast(time, event, group = g, control = 0, method = "nph")
#'
#' @references
#' Fleming, T. R., & Harrington, D. P. (1991). Counting Processes and Survival
#' Analysis. New York: John Wiley & Sons.
#'
#' @importFrom stats qnorm pnorm sd
#' @export
medsurv_fast <- function(time, event, group = NULL, control = NULL,
                         side = 2, conf.level = 0.95, conf.type = "log",
                         method = c("km", "nph"), bw = NULL,
                         presorted = FALSE) {
  method <- match.arg(method)

  n_obs <- length(time)
  if (length(event) != n_obs) {
    stop("time and event must have the same length.")
  }
  if (n_obs < 1L) {
    stop("time must contain at least one observation.")
  }
  if (anyNA(time) || anyNA(event)) {
    stop("time and event must not contain NA.")
  }
  if (!all(event %in% c(0, 1))) {
    stop("event must be coded as 0 (censored) or 1 (event).")
  }
  if (!(side %in% c(1, 2))) {
    stop("side must be 1 or 2.")
  }
  if (!(conf.type %in% c("plain", "log"))) {
    stop("conf.type must be \"plain\" or \"log\".")
  }
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be in (0, 1).")
  }

  two_group <- !is.null(group)
  control_label <- NA
  treatment_label <- NA

  if (two_group) {
    if (length(group) != n_obs) {
      stop("group must have the same length as time.")
    }
    if (anyNA(group)) {
      stop("group must not contain NA.")
    }
    levels_g <- sort(unique(group))
    if (length(levels_g) != 2L) {
      stop("group must have exactly two distinct levels for a comparison.")
    }
    if (is.null(control)) {
      stop("control must be specified (the control group level) when group is supplied.")
    }
    if (!(control %in% levels_g)) {
      stop("control must be one of the levels of group.")
    }
    control_label <- control
    treatment_label <- levels_g[levels_g != control]
    gcode <- ifelse(group == control, 0L, 1L)
    ngroup <- 2L
  } else {
    gcode <- rep(0L, n_obs)
    ngroup <- 1L
  }

  if (presorted) {
    s_time <- as.numeric(time)
    s_event <- as.integer(event)
    s_grp <- gcode
  } else {
    ord <- order(time)
    s_time <- as.numeric(time)[ord]
    s_event <- as.integer(event)[ord]
    s_grp <- gcode[ord]
  }

  bw_vec <- numeric(ngroup)
  for (g in seq_len(ngroup) - 1L) {
    evt <- s_time[s_grp == g & s_event == 1L]
    m_g <- length(evt)
    b <- if (m_g >= 2L && stats::sd(evt) > 0) {
      1.06 * stats::sd(evt) * m_g^(-1 / 5)
    } else {
      NA_real_
    }
    if (!is.finite(b) || b <= 0) {
      span <- if (m_g >= 1L) diff(range(evt)) else 0
      b <- if (span > 0) span / max(m_g, 1L)^(1 / 3) else 1
    }
    bw_vec[g + 1L] <- b
  }
  if (!is.null(bw)) {
    if (length(bw) == 1L) {
      bw_vec[] <- bw
    } else if (length(bw) == ngroup) {
      bw_vec <- bw
    } else {
      stop("bw must be a single value or have length equal to the number of groups.")
    }
    if (any(!is.finite(bw_vec)) || any(bw_vec <= 0)) {
      stop("bw must be positive and finite.")
    }
  }

  core <- medsurv_core(s_time, s_event, s_grp, ngroup, bw_vec)
  median_v <- core[, 1]
  greenwood_v <- core[, 3]
  kernel_haz_v <- core[, 4]
  nph_var_v <- core[, 5]
  local_haz_v <- core[, 6]

  if (method == "nph") {
    var_v <- ifelse(is.finite(local_haz_v) & local_haz_v > 0,
                    nph_var_v / (local_haz_v * local_haz_v), NA_real_)
  } else {
    var_v <- ifelse(is.finite(kernel_haz_v) & kernel_haz_v > 0,
                    greenwood_v / (kernel_haz_v * kernel_haz_v), NA_real_)
  }
  se_v <- sqrt(var_v)

  zc <- stats::qnorm(1 - (1 - conf.level) / 2)

  median_ci <- function(m_i, s_i) {
    if (!is.finite(m_i) || !is.finite(s_i)) {
      return(c(NA_real_, NA_real_))
    }
    if (conf.type == "log" && m_i > 0) {
      c(exp(log(m_i) - zc * s_i / m_i), exp(log(m_i) + zc * s_i / m_i))
    } else {
      c(m_i - zc * s_i, m_i + zc * s_i)
    }
  }

  n_g <- integer(ngroup)
  event_g <- integer(ngroup)
  for (g in seq_len(ngroup) - 1L) {
    n_g[g + 1L] <- sum(s_grp == g)
    event_g[g + 1L] <- sum(s_grp == g & s_event == 1L)
  }

  if (!two_group) {
    ci <- median_ci(median_v[1], se_v[1])
    out <- c(median = median_v[1], se = se_v[1], lower = ci[1], upper = ci[2])
  } else {
    ci0 <- median_ci(median_v[1], se_v[1])
    ci1 <- median_ci(median_v[2], se_v[2])
    diff_v <- median_v[2] - median_v[1]
    se_diff <- sqrt(var_v[1] + var_v[2])
    z <- diff_v / se_diff
    chisq <- z * z
    p <- if (side == 2) {
      2 * stats::pnorm(-abs(z))
    } else {
      stats::pnorm(z, lower.tail = FALSE)
    }
    diff_ci <- if (is.finite(diff_v) && is.finite(se_diff)) {
      c(diff_v - zc * se_diff, diff_v + zc * se_diff)
    } else {
      c(NA_real_, NA_real_)
    }
    out <- c(median.control = median_v[1], median.treatment = median_v[2],
             diff = diff_v,
             se.control = se_v[1], se.treatment = se_v[2], se.diff = se_diff,
             lower.control = ci0[1], upper.control = ci0[2],
             lower.treatment = ci1[1], upper.treatment = ci1[2],
             lower.diff = diff_ci[1], upper.diff = diff_ci[2],
             z = z, chisq = chisq, p = p)
  }

  attr(out, "two.group") <- two_group
  attr(out, "conf.level") <- conf.level
  attr(out, "conf.type") <- conf.type
  attr(out, "side") <- side
  attr(out, "n") <- n_g
  attr(out, "nevent") <- event_g
  attr(out, "bw") <- bw_vec
  attr(out, "control") <- control_label
  attr(out, "treatment") <- treatment_label
  attr(out, "method") <- method
  class(out) <- "medsurv_fast"
  out
}
