#' Fast Restricted Mean Survival Time (Single Group or Two-Group Comparison)
#'
#' @description
#' Computes the restricted mean survival time (RMST) up to a horizon \code{tau}
#' from the Kaplan-Meier estimator. With a single group (the default), it
#' returns the RMST with its Greenwood-type standard error and a Wald
#' confidence interval. When a \code{group} is supplied, it additionally
#' returns the two-group contrasts: the RMST difference (treatment minus
#' control) and the RMST ratio (treatment over control), each with a standard
#' error, confidence interval, and two-sided test. The C++ backend integrates
#' the survival step function in a single scan and is reused once per group, so
#' the function is suitable for simulation loops with \code{presorted = TRUE}.
#'
#' @details
#' The RMST is the area under the Kaplan-Meier curve from 0 to \code{tau}:
#'
#' RMST(tau) = integral over [0, tau] of S(u) du.
#'
#' The variance follows the Greenwood-type estimator
#'
#' Var[RMST] = sum_\{t_i <= tau\} A_i^2 d_i / (n_i (n_i - d_i)),
#'
#' where A_i = integral over [t_i, tau] of S(u) du is the area to the right of
#' event time t_i. This matches the restricted mean reported by
#' \code{\link[survival]{survfit}} and by the \code{survRM2} package.
#'
#' When \code{group} is supplied, the treatment group is the level of
#' \code{group} that is not equal to \code{control}. Writing the treatment and
#' control RMST as r1 and r0 with variances v1 and v0, the difference contrast
#' is diff = r1 - r0 with Var[diff] = v1 + v0, using the independence of the
#' two groups. The ratio contrast is formed on the log scale by the delta
#' method, Var[log(r1 / r0)] = v1 / r1^2 + v0 / r0^2, with the confidence
#' interval exponentiated back to the ratio scale. These match the unadjusted
#' contrasts reported by \code{survRM2}.
#'
#' When \code{presorted = TRUE}, the input vectors are assumed to be sorted in
#' ascending order of \code{time}; splitting into groups preserves the
#' ascending order within each group, so no re-sorting is performed. When
#' \code{presorted = FALSE} (default), sorting is handled internally. In
#' simulation loops where the data are generated in sorted order,
#' \code{presorted = TRUE} avoids one O(n log n) pass.
#'
#' @param time A numeric vector of event or censoring times.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group An optional vector of group labels aligned with \code{time}.
#'   If \code{NULL} (default), a single-group RMST is computed. If supplied,
#'   the two-group contrasts are returned and \code{control} must be given.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group. Required when \code{group} is supplied.
#' @param side 1 for a one-sided test in the direction of treatment benefit
#'   (longer restricted mean survival in the treatment group) or 2 for a
#'   two-sided test (default 2). Applies to the two-group contrasts; the
#'   confidence intervals are always two-sided at \code{conf.level}.
#' @param conf.level A single numeric value in (0, 1) specifying the confidence
#'   level. Defaults to 0.95.
#' @param tau A single positive numeric value specifying the restriction
#'   horizon.
#' @param presorted A logical value. If \code{TRUE}, the inputs are assumed to
#'   be sorted in ascending order of \code{time}. If \code{FALSE} (default),
#'   sorting is handled internally.
#'
#' @return An object of class \code{"rmst_fast"}, a named numeric vector. In
#'   single-group mode it has length 4 with elements \code{rmst},
#'   \code{std.err}, \code{lower}, and \code{upper}. In two-group mode it
#'   contains the per-group RMST (\code{rmst.ctrl}, \code{rmst.trt}), the
#'   difference contrast (\code{diff}, \code{se.diff}, \code{diff.lower},
#'   \code{diff.upper}, \code{z.diff}, \code{p.diff}), and the ratio contrast
#'   (\code{ratio}, \code{ratio.lower}, \code{ratio.upper}, \code{z.ratio},
#'   \code{p.ratio}). The restriction horizon, confidence level, and test side
#'   are stored as attributes \code{tau}, \code{conf.level}, and \code{side};
#'   in two-group mode the \code{control} label is also stored. Returns
#'   \code{NA_real_} values (still with class \code{"rmst_fast"}) when \code{n}
#'   is zero in single-group mode.
#'
#' @examples
#' set.seed(42)
#' t_raw <- rexp(100, rate = 1 / 10)
#' e_raw <- rbinom(100, 1, 0.7)
#'
#' # Single-group RMST
#' rmst_fast(t_raw, e_raw, tau = 10)
#'
#' # Single-group, pre-sorted (sort once, reuse in a loop)
#' ord <- order(t_raw)
#' rmst_fast(t_raw[ord], e_raw[ord], tau = 10, presorted = TRUE)
#'
#' # Two-group comparison (difference and ratio)
#' set.seed(7)
#' n <- 200
#' time  <- c(rexp(n, 0.10), rexp(n, 0.07))
#' event <- rbinom(2 * n, 1, 0.8)
#' group <- rep(0:1, each = n)
#' rmst_fast(time, event, group = group, control = 0, tau = 10)
#'
#' \donttest{
#' # Validation against survRM2
#' if (requireNamespace("survRM2", quietly = TRUE)) {
#'   survRM2::rmst2(time, event, group, tau = 10)$unadjusted.result
#' }
#' }
#'
#' @references
#' Royston, P., & Parmar, M. K. B. (2013). Restricted mean survival time: an
#' alternative to the hazard ratio for the design and analysis of randomized
#' trials with a time-to-event outcome. \emph{BMC Medical Research
#' Methodology}, \emph{13}, 152.
#'
#' Uno, H., Claggett, B., Tian, L., et al. (2014). Moving beyond the hazard
#' ratio in quantifying the between-group difference in survival analysis.
#' \emph{Journal of Clinical Oncology}, \emph{32}(22), 2380-2385.
#'
#' @seealso
#' \code{\link{survfit_fast}} for the Kaplan-Meier estimate at a single time
#' point.
#'
#' @importFrom stats qnorm pnorm
#' @export
rmst_fast <- function(time, event, group = NULL, control = NULL, side = 2,
                      conf.level = 0.95, tau, presorted = FALSE) {
  # Input validation
  n <- length(time)
  if (length(event) != n) {
    stop("'time' and 'event' must have the same length")
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

  # ---- Single-group mode -------------------------------------------------
  if (is.null(group)) {
    na_out <- c(rmst = NA_real_, std.err = NA_real_,
                lower = NA_real_, upper = NA_real_)
    wrap1  <- function(v) {
      structure(v, tau = tau, conf.level = conf.level, side = side,
                class = "rmst_fast")
    }
    if (n == 0L) return(wrap1(na_out))

    if (!presorted) {
      ord   <- order(time)
      time  <- time[ord]
      event <- event[ord]
    }

    res      <- rmst_core(time, event, tau)
    rmst     <- res[1L]
    var_rmst <- res[2L]
    if (is.na(rmst)) return(wrap1(na_out))

    std.err <- sqrt(var_rmst)
    return(wrap1(c(rmst    = rmst,
                   std.err = std.err,
                   lower   = rmst - z * std.err,
                   upper   = rmst + z * std.err)))
  }

  # ---- Two-group mode ----------------------------------------------------
  if (length(group) != n) {
    stop("'group' must have the same length as 'time'")
  }
  if (is.null(control)) {
    stop("'control' must be supplied when 'group' is given")
  }

  # Treatment indicator: 1 = treatment, 0 = control
  if (is.factor(group)) group <- as.character(group)
  j <- as.integer(group != control)

  # Sort pooled data by time when not presorted (group split preserves order)
  if (!presorted) {
    ord   <- order(time)
    time  <- time[ord]
    event <- event[ord]
    j     <- j[ord]
  }

  is1  <- j == 1L
  res0 <- rmst_core(time[!is1], event[!is1], tau)
  res1 <- rmst_core(time[is1],  event[is1],  tau)

  r0 <- res0[1L]; v0 <- res0[2L]
  r1 <- res1[1L]; v1 <- res1[2L]

  # Difference contrast (treatment minus control)
  est_diff <- r1 - r0
  se_diff  <- sqrt(v1 + v0)
  z_diff   <- if (is.finite(se_diff) && se_diff > 0) est_diff / se_diff else NA_real_
  # side = 1 tests the treatment-benefit direction (longer RMST in treatment,
  # i.e. a positive contrast, upper tail).
  p_diff   <- if (is.na(z_diff)) {
    NA_real_
  } else if (side == 1L) {
    pnorm(-z_diff)
  } else {
    2 * pnorm(-abs(z_diff))
  }

  # Ratio contrast on the log scale (delta method)
  if (is.finite(r0) && is.finite(r1) && r0 > 0 && r1 > 0) {
    ratio    <- r1 / r0
    se_lr    <- sqrt(v1 / (r1 * r1) + v0 / (r0 * r0))
    log_r    <- log(ratio)
    z_ratio  <- if (se_lr > 0) log_r / se_lr else NA_real_
    p_ratio  <- if (is.na(z_ratio)) {
      NA_real_
    } else if (side == 1L) {
      pnorm(-z_ratio)
    } else {
      2 * pnorm(-abs(z_ratio))
    }
    ratio_lo <- exp(log_r - z * se_lr)
    ratio_hi <- exp(log_r + z * se_lr)
  } else {
    ratio    <- NA_real_; ratio_lo <- NA_real_; ratio_hi <- NA_real_
    z_ratio  <- NA_real_; p_ratio  <- NA_real_
  }

  out <- c(
    rmst.ctrl   = r0,
    rmst.trt    = r1,
    diff        = est_diff,
    se.diff     = se_diff,
    diff.lower  = est_diff - z * se_diff,
    diff.upper  = est_diff + z * se_diff,
    z.diff      = z_diff,
    p.diff      = p_diff,
    ratio       = ratio,
    ratio.lower = ratio_lo,
    ratio.upper = ratio_hi,
    z.ratio     = z_ratio,
    p.ratio     = p_ratio
  )

  structure(out, tau = tau, conf.level = conf.level, side = side,
            control = control, class = "rmst_fast")
}
