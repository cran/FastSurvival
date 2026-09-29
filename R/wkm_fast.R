#' Fast weighted Kaplan-Meier (Pepe-Fleming) two-sample test
#'
#' Computes the weighted Kaplan-Meier (WKM) statistic of Pepe and Fleming for
#' comparing two survival curves. The statistic is the weighted integral of the
#' difference between the Kaplan-Meier estimates of the treatment and control
#' groups over the observed range, standardized to a Wald z statistic. Unlike
#' weighted log-rank tests, this test targets the integrated difference in
#' survival and is sensitive to differences even when the hazard functions
#' cross. The single scan over the sorted data is performed in C++ for use
#' inside simulation loops.
#'
#' The default weight is the Pepe-Fleming combined censoring weight
#' w(t) = n G1(t) G2(t) / (n1 G1(t) + n2 G2(t)), where G1 and G2 are the
#' Kaplan-Meier estimates of the censoring survival functions of the two
#' groups. This weight gives less weight to regions with heavy censoring and
#' stabilizes the variance in the
#' tail. The choice \code{weight = "sqrtPF"} uses its square root, and
#' \code{weight = "constant"} uses a weight of 1, in which case the numerator
#' reduces to the difference in restricted mean survival time over the observed
#' range. With \code{weight = "PF"} the result reproduces
#' \code{nphsim::wkm.Stat} for data without tied times.
#'
#' The weighted difference is computed as treatment minus control, so a
#' positive value and a positive z indicate longer survival under treatment.
#'
#' @param time Numeric vector of event or censoring times.
#' @param event Integer vector, 1 for an event and 0 for censoring.
#' @param group Grouping vector with exactly two distinct levels.
#' @param control The level of \code{group} that denotes the control group.
#' @param side Either 2 for a two-sided test or 1 for a one-sided test of
#'   treatment superiority (weighted difference greater than 0).
#' @param conf.level Confidence level for the interval of the weighted
#'   difference.
#' @param weight Weight function, one of \code{"PF"} (Pepe-Fleming combined
#'   censoring weight), \code{"sqrtPF"} (its square root) or \code{"constant"}
#'   (weight 1).
#' @param presorted Logical; set to \code{TRUE} when \code{time}, \code{event}
#'   and \code{group} are already sorted in ascending order of \code{time}, to
#'   skip the internal sort. Defaults to \code{FALSE}.
#'
#' @return A named numeric vector of class \code{"wkm_fast"} with the weighted
#'   integrated difference, its standard error and confidence limits, and the
#'   Wald statistics for the test.
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' g <- rep(0:1, each = n / 2)
#' tt <- c(rexp(n / 2, log(2) / 12), rexp(n / 2, log(2) / 16))
#' cc <- rexp(n, rate = 0.02)
#' time <- pmin(tt, cc)
#' event <- as.integer(tt <= cc)
#' wkm_fast(time, event, group = g, control = 0)
#'
#' @references
#' Pepe, M. S., & Fleming, T. R. (1989). Weighted Kaplan-Meier statistics: a
#' class of distance tests for censored survival data. Biometrics, 45(2),
#' 497-507.
#'
#' Pepe, M. S., & Fleming, T. R. (1991). Weighted Kaplan-Meier statistics:
#' large sample and optimality considerations. Journal of the Royal Statistical
#' Society. Series B (Methodological), 53(2), 341-352.
#'
#' @importFrom stats qnorm pnorm
#' @export
wkm_fast <- function(time, event, group, control = NULL,
                     side = 2, conf.level = 0.95,
                     weight = c("PF", "sqrtPF", "constant"),
                     presorted = FALSE) {
  weight <- match.arg(weight)

  n_obs <- length(time)
  if (length(event) != n_obs) {
    stop("time and event must have the same length.")
  }
  if (length(group) != n_obs) {
    stop("group must have the same length as time.")
  }
  if (n_obs < 1L) {
    stop("time must contain at least one observation.")
  }
  if (anyNA(time) || anyNA(event) || anyNA(group)) {
    stop("time, event and group must not contain NA.")
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

  levels_g <- sort(unique(group))
  if (length(levels_g) != 2L) {
    stop("group must have exactly two distinct levels.")
  }
  if (is.null(control)) {
    stop("control must be specified (the control group level).")
  }
  if (!(control %in% levels_g)) {
    stop("control must be one of the levels of group.")
  }
  treatment_label <- levels_g[levels_g != control]
  gcode <- ifelse(group == control, 0L, 1L)

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

  weight_type <- switch(weight, PF = 0L, sqrtPF = 1L, constant = 2L)

  core <- wkm_core(s_time, s_event, s_grp, weight_type)
  num_raw <- core[1]
  variance <- core[2]
  n1 <- core[3]
  n2 <- core[4]
  n_total <- n1 + n2

  scale <- sqrt(n1 * n2 / n_total)
  if (is.finite(variance) && variance > 0) {
    se_wdiff <- sqrt(variance) / scale
    z <- num_raw / se_wdiff
  } else {
    se_wdiff <- NA_real_
    z <- NA_real_
  }
  chisq <- z * z
  p <- if (side == 2) {
    2 * stats::pnorm(-abs(z))
  } else {
    stats::pnorm(z, lower.tail = FALSE)
  }

  zc <- stats::qnorm(1 - (1 - conf.level) / 2)
  lower <- if (is.finite(num_raw) && is.finite(se_wdiff)) num_raw - zc * se_wdiff else NA_real_
  upper <- if (is.finite(num_raw) && is.finite(se_wdiff)) num_raw + zc * se_wdiff else NA_real_

  n_g <- c(n2, n1)
  event_g <- c(sum(s_grp == 0 & s_event == 1L), sum(s_grp == 1 & s_event == 1L))

  out <- c(wdiff = num_raw, se = se_wdiff, lower = lower, upper = upper,
           z = z, chisq = chisq, p = p)

  attr(out, "conf.level") <- conf.level
  attr(out, "side") <- side
  attr(out, "weight") <- weight
  attr(out, "n") <- n_g
  attr(out, "nevent") <- event_g
  attr(out, "control") <- control
  attr(out, "treatment") <- treatment_label
  class(out) <- "wkm_fast"
  out
}
