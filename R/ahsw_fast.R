#' Fast Average Hazard with Survival Weight (Two-Group Comparison)
#'
#' @description
#' Computes the average hazard with survival weight (AHSW) of Uno and Horiguchi
#' for two groups and the between-group contrasts. The average hazard on the
#' window from 0 to \code{tau} is the ratio of the cumulative event probability
#' at \code{tau} to the restricted mean survival time at \code{tau}, both based
#' on the Kaplan-Meier estimate. The function returns the per-group average
#' hazard, the ratio of average hazards (RAH, treatment over control) on the log
#' scale, and the difference of average hazards (DAH, treatment minus control)
#' on the identity scale, each with a confidence interval and a two-sided test.
#' The C++ backend walks the pooled sorted data once per group, so the function
#' is suitable for simulation loops with \code{presorted = TRUE}.
#'
#' @details
#' The average hazard with survival weight is
#'
#' AH(tau) = (1 - S(tau)) / integral over [0, tau] of S(u) du,
#'
#' the ratio of the cumulative event probability at \code{tau} to the restricted
#' mean survival time at \code{tau}. It can be read as a general censoring-free
#' incidence rate on the window from 0 to \code{tau} and stays interpretable
#' under non-proportional hazards. This is a different quantity from the average
#' hazard ratio of Kalbfleisch, which averages the time-varying ratio of hazards
#' rather than forming a single average hazard per group and then contrasting.
#'
#' Writing the treatment and control average hazards as a1 and a0, the ratio
#' contrast is RAH = a1 / a0, formed on the log scale with variance
#' v_Q1 / n1 + v_Q0 / n0, and the difference contrast is DAH = a1 - a0, with
#' variance v_U1 / n1 + v_U0 / n0, using the independence of the two groups.
#' The per-group variance terms v_Q (log scale) and v_U (identity scale) follow
#' the asymptotic variance of Uno and Horiguchi, computed from the Nelson-Aalen
#' increments, the running restricted mean survival time and the at-risk
#' fraction. The confidence interval for RAH is exponentiated from the log
#' scale, and the two-sided p-values are based on the normal approximation.
#'
#' When \code{presorted = TRUE}, the inputs are assumed to be sorted in
#' ascending order of \code{time}, so the internal \code{order()} call is
#' skipped. Splitting into groups preserves the ascending order within each
#' group.
#'
#' @param time A numeric vector of follow-up times for all subjects.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param side 1 for a one-sided test in the direction of treatment benefit
#'   (lower average hazard in the treatment group, i.e. a ratio below 1 and a
#'   negative difference) or 2 for a two-sided test (default 2). The reported
#'   p-values follow this choice; the confidence intervals are always two-sided
#'   at \code{conf.level}.
#' @param conf.level A single numeric value in (0, 1) specifying the confidence
#'   level. Defaults to 0.95.
#' @param tau A single positive numeric value, the truncation time point for the
#'   average hazard. Both groups must have positive Kaplan-Meier survival at
#'   \code{tau}.
#' @param presorted A logical value. If \code{TRUE}, \code{time}, \code{event},
#'   and \code{group} are assumed to be sorted in ascending order of
#'   \code{time}, and the internal \code{order()} call is skipped. If
#'   \code{FALSE} (default), sorting is handled internally.
#'
#' @return An object of class \code{"ahsw_fast"}, a named numeric vector
#'   containing the per-group average hazards (\code{ah.ctrl}, \code{ah.trt}),
#'   the ratio contrast (\code{rah}, \code{rah.lower}, \code{rah.upper},
#'   \code{p.rah}), and the difference contrast (\code{dah}, \code{dah.lower},
#'   \code{dah.upper}, \code{p.dah}). The truncation time and confidence level
#'   are stored as attributes \code{tau} and \code{conf.level}, and the
#'   \code{control} label is also stored. Returns \code{NA} values (still with
#'   class \code{"ahsw_fast"}) when either group has zero survival at \code{tau}
#'   or a non-finite variance.
#'
#' @examples
#' library(survival)
#'
#' # Average hazard contrasts on the ovarian data
#' ahsw_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1, tau = 600)
#'
#' # presorted = TRUE: sort once outside, reuse inside a loop
#' ord <- order(ovarian$futime)
#' ahsw_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
#'           control = 1, tau = 600, presorted = TRUE)
#'
#' \donttest{
#' # Validation against survAH
#' if (requireNamespace("survAH", quietly = TRUE)) {
#'   arm <- as.numeric(ovarian$rx == 2)
#'   survAH::ah2(time = ovarian$futime, status = ovarian$fustat,
#'               arm = arm, tau = 600)
#' }
#' }
#'
#' @references
#' Uno, H., & Horiguchi, M. (2023). Ratio and difference of average hazard with
#' survival weight: new measures to quantify survival benefit of new therapy.
#' \emph{Statistics in Medicine}, \emph{42}(7), 936-952.
#'
#' @seealso
#' \code{\link{rmst_fast}} for the restricted mean survival time, and
#' \code{\link{survfit_fast}} for the Kaplan-Meier estimate at a time point.
#'
#' @importFrom stats qnorm pnorm
#' @export
ahsw_fast <- function(time, event, group, control, side = 2,
                      conf.level = 0.95, tau, presorted = FALSE) {
  # Input validation
  n_total <- length(time)
  if (length(event) != n_total || length(group) != n_total) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("'tau' must be a single positive value")
  }
  if (conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be in (0, 1)")
  }

  z <- qnorm(1 - (1 - conf.level) / 2)

  # Treatment indicator: 1 = treatment, 0 = control
  if (is.factor(group)) group <- as.character(group)
  j <- as.integer(group != control)

  # Sort pooled data by time when not presorted (group split preserves order)
  if (!presorted) {
    ord   <- order(time)
    time  <- time[ord]
    event <- as.integer(event[ord])
    j     <- j[ord]
  } else {
    event <- as.integer(event)
  }

  is1  <- j == 1L
  n0   <- sum(!is1)
  n1   <- sum(is1)

  res0 <- ahsw_core(time[!is1], event[!is1], tau)
  res1 <- ahsw_core(time[is1],  event[is1],  tau)

  a0 <- res0[3L]; vQ0 <- res0[4L]; vU0 <- res0[5L]
  a1 <- res1[3L]; vQ1 <- res1[4L]; vU1 <- res1[5L]

  out <- c(
    ah.ctrl   = a0,
    ah.trt    = a1,
    rah       = NA_real_,
    rah.lower = NA_real_,
    rah.upper = NA_real_,
    p.rah     = NA_real_,
    dah       = NA_real_,
    dah.lower = NA_real_,
    dah.upper = NA_real_,
    p.dah     = NA_real_
  )

  wrap <- function(v) {
    structure(v, tau = tau, conf.level = conf.level, side = side,
              control = control, class = "ahsw_fast")
  }

  # Bail out if either group has zero survival at tau or a non-finite variance
  if (!is.finite(a0) || !is.finite(a1) || a0 <= 0 || a1 <= 0 ||
      !is.finite(vQ0) || !is.finite(vQ1) ||
      !is.finite(vU0) || !is.finite(vU1)) {
    return(wrap(out))
  }

  # Ratio of average hazards on the log scale. Treatment benefit is a ratio
  # below 1 (log_rah < 0); the one-sided test (side = 1) uses the lower tail.
  log_rah <- log(a1 / a0)
  se_rah  <- sqrt(vQ1 / n1 + vQ0 / n0)
  rah     <- exp(log_rah)
  p_rah   <- if (se_rah <= 0) {
    NA_real_
  } else if (side == 1L) {
    pnorm(log_rah / se_rah)
  } else {
    2 * pnorm(-abs(log_rah) / se_rah)
  }

  # Difference of average hazards on the identity scale. Treatment benefit is a
  # negative difference; the one-sided test uses the lower tail.
  dah    <- a1 - a0
  se_dah <- sqrt(vU1 / n1 + vU0 / n0)
  p_dah  <- if (se_dah <= 0) {
    NA_real_
  } else if (side == 1L) {
    pnorm(dah / se_dah)
  } else {
    2 * pnorm(-abs(dah) / se_dah)
  }

  out["rah"]       <- rah
  out["rah.lower"] <- exp(log_rah - z * se_rah)
  out["rah.upper"] <- exp(log_rah + z * se_rah)
  out["p.rah"]     <- p_rah
  out["dah"]       <- dah
  out["dah.lower"] <- dah - z * se_dah
  out["dah.upper"] <- dah + z * se_dah
  out["p.dah"]     <- p_dah

  wrap(out)
}
