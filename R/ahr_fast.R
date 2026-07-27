#' Fast Kalbfleisch-Prentice average hazard ratio for two groups
#'
#' Estimates the average hazard ratio of Kalbfleisch and Prentice (1981) between
#' two groups over the time interval from 0 to \code{tau}, based on the
#' Kaplan-Meier estimator of each group's survival function. This is the
#' Kaplan-Meier (unweighted) special case of the estimator implemented in the
#' archived \code{AHR} package, restricted to two groups and recoded in C++ for
#' use inside simulation loops.
#'
#' The estimator works with the group shares of the total hazard. Writing
#' \code{S1} and \code{S2} for the two survival functions, the reference-group
#' share is \code{theta1 = -integral(S2 dS1) / (1 - S1(tau) S2(tau))} over
#' \code{[0, tau]}, the comparison-group share is \code{theta2 = 1 - theta1},
#' and the average hazard ratio is \code{ahr = theta2 / theta1}. Under
#' proportional hazards with hazard ratio \code{psi}, \code{ahr} estimates
#' \code{psi}. A value above 1 indicates higher hazard (worse survival) in the
#' comparison group. The variance of \code{theta1} is the direct
#' Greenwood-based estimator. The primary test is on the theta (group-share)
#' scale, as in Kalbfleisch and Prentice (1981) and Dormuth et al. (2024,
#' eq. 5): the comparison-group share is compared with its null value (0.5 when
#' \code{null.ahr = 1}). An equivalent test and a confidence interval on the
#' \code{log(ahr)} scale are also reported.
#'
#' This is distinct from \code{\link{ahsw_fast}}, which estimates the
#' Uno-Horiguchi average hazard with survival weight.
#'
#' @param time vector of right-censored event times
#' @param event 0/1 (or logical) event indicators, 1 for an event
#' @param group vector with exactly two distinct values identifying the groups.
#' @param control the value of \code{group} that denotes the reference (control)
#'   group. The average hazard ratio is reported for the other group relative to
#'   it.
#' @param side 1 for a one-sided test in the direction of treatment benefit
#'   (average hazard ratio below \code{null.ahr}) or 2 for a two-sided test
#'   (default 2).
#' @param conf.level confidence level for the confidence interval (default 0.95)
#' @param tau upper limit of the interval over which the average hazard ratio is
#'   computed. If \code{NULL} (default) the largest time observed in both groups
#'   is used.
#' @param null.ahr value of the average hazard ratio under the null hypothesis
#'   used for the Z statistic and p-value (default 1)
#' @param presorted if \code{TRUE}, assume \code{time} is already sorted in
#'   ascending order so that each group's observations are also ascending; this
#'   skips the internal sort (default \code{FALSE})
#' @return An object of class \code{"ahr_fast"}, a list with elements
#'   \code{ahr} (the average hazard ratio, comparison vs reference),
#'   \code{log.ahr}, \code{se.loghr}, \code{lower}, \code{upper},
#'   \code{conf.level}, \code{z} and \code{p.value} (the primary test on the
#'   theta / group-share scale, as in Dormuth et al. 2024 eq. 5),
#'   \code{z.loghr} and \code{p.value.loghr} (the equivalent test on the
#'   \code{log(ahr)} scale), \code{se.theta} (standard error of the tested
#'   comparison-group share), \code{null.share}, \code{null.ahr}, \code{theta}
#'   (the two group shares), \code{var.theta1}, \code{var.theta2},
#'   \code{side}, \code{tau}, \code{n} (the two group sizes) and \code{groups}.
#' @references
#' Kalbfleisch, J. D., & Prentice, R. L. (1981). Estimation of the average
#' hazard ratio. \emph{Biometrika}, \emph{68}(1), 105-112.
#'
#' Dormuth, I., Pauly, M., Rauch, G., & Herrmann, C. (2024). Sample size
#' calculation under nonproportional hazards using average hazard ratios.
#' \emph{Biometrical Journal}, \emph{66}(6), e202300271.
#' @seealso \code{\link{ahsw_fast}}, \code{\link{coxph_fast}}
#' @export
#' @examples
#' set.seed(1)
#' n <- 200
#' time1 <- rexp(n, 0.1)
#' time2 <- rexp(n, 0.18)
#' cens <- rexp(2 * n, 0.05)
#' obs <- pmin(c(time1, time2), cens)
#' event <- as.integer(c(time1, time2) <= cens)
#' group <- rep(c(0, 1), each = n)
#' ahr_fast(obs, event, group, control = 0, tau = 8)
ahr_fast <- function(time, event, group, control, side = 2,
                     conf.level = 0.95, tau = NULL, null.ahr = 1,
                     presorted = FALSE) {

  # Input validation
  if (!is.numeric(time)) stop("'time' must be numeric.")
  event <- as.integer(event)
  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided).")
  }
  if (length(time) != length(event) || length(time) != length(group)) {
    stop("'time', 'event' and 'group' must have the same length.")
  }
  if (anyNA(time) || anyNA(event) || anyNA(group)) {
    stop("'time', 'event' and 'group' must not contain missing values.")
  }
  if (any(!(event %in% c(0L, 1L)))) stop("'event' must be 0 or 1.")
  if (any(time < 0)) stop("'time' must be non-negative.")
  if (!is.numeric(null.ahr) || length(null.ahr) != 1 || null.ahr <= 0) {
    stop("'null.ahr' must be a single positive number.")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1 ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number in (0, 1).")
  }

  lev <- sort(unique(group))
  if (length(lev) != 2) stop("'group' must have exactly two distinct values.")
  if (missing(control) || length(control) != 1 ||
      !(as.character(control) %in% as.character(lev))) {
    stop("'control' must be one of the two values in 'group'.")
  }
  ref <- lev[as.character(lev) == as.character(control)]
  cmp <- lev[as.character(lev) != as.character(control)]

  sel1 <- group == ref
  sel2 <- group == cmp
  time1 <- time[sel1]
  event1 <- event[sel1]
  time2 <- time[sel2]
  event2 <- event[sel2]
  n1 <- length(time1)
  n2 <- length(time2)
  if (n1 == 0 || n2 == 0) stop("Both groups must be non-empty.")

  if (!presorted) {
    ord1 <- order(time1)
    ord2 <- order(time2)
    time1 <- time1[ord1]
    event1 <- event1[ord1]
    time2 <- time2[ord2]
    event2 <- event2[ord2]
  }

  if (is.null(tau)) {
    tau <- min(max(time1), max(time2))
  } else {
    if (!is.numeric(tau) || length(tau) != 1 || tau <= 0) {
      stop("'tau' must be a single positive number or NULL.")
    }
  }

  # Common evaluation grid: 0, distinct event times up to tau, and tau
  events <- c(time1[event1 == 1], time2[event2 == 1])
  events <- events[events <= tau]
  eval_times <- sort(unique(c(0, events, tau)))

  core <- ahr_core(time1, event1, time2, event2, eval_times)

  theta1 <- core[["theta1"]]
  theta2 <- core[["theta2"]]
  ahr <- core[["ahr"]]
  var.theta1 <- core[["var.theta1"]]
  var.theta2 <- core[["var.theta2"]]
  valid <- core[["valid"]] == 1

  if (valid) {
    # Primary test on the theta (group-share) scale, as in Kalbfleisch and
    # Prentice (1981) and Dormuth et al. (2024, eq. 5). Under H0 that
    # AHR = null.ahr, the comparison-group share equals null.ahr / (1 +
    # null.ahr) (which is 0.5 when null.ahr = 1). Using the comparison share
    # (rather than the reference share) aligns the sign of z with log(ahr).
    # The comparison share is tested with its own variance (var.theta2), which
    # matches the per-component test in the AHR package (the direct estimator
    # does not force var.theta1 == var.theta2 in finite samples).
    se.theta <- sqrt(var.theta2)
    null.share <- null.ahr / (1 + null.ahr)
    z <- (theta2 - null.share) / se.theta
    # Treatment benefit corresponds to an average hazard ratio below null.ahr,
    # i.e. a comparison-group share below null.share, hence a negative z. A
    # one-sided test (side = 1) uses the lower tail.
    p.value <- if (side == 1L) {
      stats::pnorm(z)
    } else {
      2 * stats::pnorm(-abs(z))
    }

    # Average hazard ratio effect: log-scale confidence interval and the
    # equivalent test on the log(ahr) scale (using var.theta1, as the AHR
    # package's log-hazard-ratio inference does)
    log.ahr <- log(ahr)
    se.loghr <- sqrt(var.theta1 / (theta1 * theta2)^2)
    z.loghr <- (log.ahr - log(null.ahr)) / se.loghr
    p.value.loghr <- if (side == 1L) {
      stats::pnorm(z.loghr)
    } else {
      2 * stats::pnorm(-abs(z.loghr))
    }
    zq <- stats::qnorm(1 - (1 - conf.level) / 2)
    lower <- exp(log.ahr - zq * se.loghr)
    upper <- exp(log.ahr + zq * se.loghr)
  } else {
    se.theta <- null.share <- z <- p.value <- log.ahr <- se.loghr <-
      z.loghr <- p.value.loghr <- lower <- upper <- NA_real_
    warning("Average hazard ratio is undefined for these data (no events in ",
            "[0, tau] or a degenerate group share); returning NA.")
  }

  theta <- c(theta1, theta2)
  names(theta) <- as.character(c(ref, cmp))
  n <- c(n1, n2)
  names(n) <- as.character(c(ref, cmp))

  out <- list(ahr = ahr, log.ahr = log.ahr, se.loghr = se.loghr,
              lower = lower, upper = upper, conf.level = conf.level,
              z = z, p.value = p.value, se.theta = se.theta,
              null.share = null.share, z.loghr = z.loghr,
              p.value.loghr = p.value.loghr, null.ahr = null.ahr,
              theta = theta, var.theta1 = var.theta1, var.theta2 = var.theta2,
              side = side, tau = tau, n = n, groups = c(ref, cmp))
  class(out) <- "ahr_fast"
  out
}
