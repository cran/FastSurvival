#' Fast window mean survival time and between-group difference
#'
#' Estimates the window mean survival time (WMST) for a single group, or the
#' difference in WMST between a treatment group and a control group, together
#' with standard errors, confidence intervals and a Wald test for the
#' difference. The WMST is the area under the Kaplan-Meier curve between the
#' window limits \code{tau1} and \code{tau2}. With \code{tau1 = 0} it is the
#' restricted mean survival time at \code{tau2}, and the result then matches
#' \code{survRM2::rmst2}. The single scan over the sorted data is performed in
#' C++ for use inside simulation loops.
#'
#' The variance of the WMST estimate uses the counting process (Greenwood) form
#' in which each event time at or before \code{tau2} contributes the squared
#' remaining window area. Event times at or before \code{tau1} contribute the
#' full window area, since the hazard before \code{tau1} scales the survival
#' curve throughout the window.
#'
#' The difference is computed as treatment minus control, so a positive
#' difference indicates a longer mean survival time within the window under
#' treatment.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Integer vector, 1 for an event and 0 for censoring.
#' @param group Optional grouping vector with exactly two distinct levels for a
#'   two-group comparison. If omitted, a single-group WMST is returned.
#' @param control The level of \code{group} that denotes the control group.
#'   Required when \code{group} is supplied.
#' @param tau1 Lower limit of the window. Defaults to 0.
#' @param tau2 Upper limit of the window. If omitted, the largest time common
#'   to both groups is used (the minimum over groups of the maximum observed
#'   time), or the maximum observed time for a single group. A supplied value
#'   larger than that gives a warning, because the Kaplan-Meier curve is not
#'   estimated beyond it.
#' @param side Either 2 for a two-sided test or 1 for a one-sided test of
#'   treatment superiority (difference greater than 0).
#' @param conf.level Confidence level for the intervals.
#' @param presorted Logical; set to \code{TRUE} when \code{time}, \code{event}
#'   and \code{group} are already sorted in ascending order of \code{time}, to
#'   skip the internal sort. Defaults to \code{FALSE}.
#'
#' @return A named numeric vector of class \code{"wmst_fast"}. For a single
#'   group the elements are the WMST, its standard error and confidence limits.
#'   For two groups the elements are the control and treatment WMST values,
#'   their difference, the standard errors, the confidence limits and the Wald
#'   statistics for the difference.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' g <- rep(0:1, each = n / 2)
#' tt <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.07))
#' cc <- rexp(n, rate = 0.02)
#' time <- pmin(tt, cc)
#' event <- as.integer(tt <= cc)
#' wmst_fast(time, event, group = g, control = 0, tau1 = 2, tau2 = 12)
#'
#' @references
#' Paukner, M., & Chappell, R. (2021). Window mean survival time. Statistics in
#' Medicine, 40(25), 5521-5533.
#'
#' @importFrom stats qnorm pnorm
#' @export
wmst_fast <- function(time, event, group = NULL, control = NULL,
                      tau1 = 0, tau2 = NULL, side = 2, conf.level = 0.95,
                      presorted = FALSE) {
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
  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be in (0, 1).")
  }
  if (tau1 < 0) {
    stop("tau1 must be non-negative.")
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

  if (is.null(tau2)) {
    if (two_group) {
      tau2 <- min(max(time[gcode == 0]), max(time[gcode == 1]))
    } else {
      tau2 <- max(time)
    }
  } else {
    check_tau_follow_up(time, if (two_group) gcode else NULL, tau2,
                        arg = "tau2")
  }
  if (tau2 <= tau1) {
    stop("tau2 must be greater than tau1.")
  }

  if (presorted) {
    s_time <- as.numeric(time)
    s_event <- as.integer(event)
    s_grp <- as.integer(gcode)
  } else {
    ord <- order(time)
    s_time <- as.numeric(time)[ord]
    s_event <- as.integer(event)[ord]
    s_grp <- as.integer(gcode)[ord]
  }

  core <- wmst_core(s_time, s_event, s_grp, ngroup, tau1, tau2)
  wmst_v <- core[, 1]
  var_v <- core[, 2]
  se_v <- sqrt(var_v)

  zc <- stats::qnorm(1 - (1 - conf.level) / 2)

  n_g <- integer(ngroup)
  event_g <- integer(ngroup)
  for (g in seq_len(ngroup) - 1L) {
    n_g[g + 1L] <- sum(s_grp == g)
    event_g[g + 1L] <- sum(s_grp == g & s_event == 1L & s_time <= tau2)
  }

  if (!two_group) {
    lower <- wmst_v[1] - zc * se_v[1]
    upper <- wmst_v[1] + zc * se_v[1]
    out <- c(wmst = wmst_v[1], se = se_v[1], lower = lower, upper = upper)
  } else {
    ci0 <- c(wmst_v[1] - zc * se_v[1], wmst_v[1] + zc * se_v[1])
    ci1 <- c(wmst_v[2] - zc * se_v[2], wmst_v[2] + zc * se_v[2])
    diff_v <- wmst_v[2] - wmst_v[1]
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
    out <- c(wmst.control = wmst_v[1], wmst.treatment = wmst_v[2],
             diff = diff_v,
             se.control = se_v[1], se.treatment = se_v[2], se.diff = se_diff,
             lower.control = ci0[1], upper.control = ci0[2],
             lower.treatment = ci1[1], upper.treatment = ci1[2],
             lower.diff = diff_ci[1], upper.diff = diff_ci[2],
             z = z, chisq = chisq, p = p)
  }

  attr(out, "two.group") <- two_group
  attr(out, "tau1") <- tau1
  attr(out, "tau2") <- tau2
  attr(out, "conf.level") <- conf.level
  attr(out, "side") <- side
  attr(out, "n") <- n_g
  attr(out, "nevent") <- event_g
  attr(out, "control") <- control_label
  attr(out, "treatment") <- treatment_label
  attr(out, "method") <- "Kaplan-Meier"
  class(out) <- "wmst_fast"
  out
}
