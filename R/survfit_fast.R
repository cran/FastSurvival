#' Fast Kaplan-Meier Survival Probability at a Specified Time Point
#'
#' @description
#' Computes the Kaplan-Meier survival probability at a specified time point,
#' together with a standard error and confidence interval based on Greenwood's
#' variance formula. The C++ backend performs binary search for the evaluation
#' cutoff and accumulates the Kaplan-Meier product and Greenwood sum in a
#' single scan over event positions only, without constructing intermediate
#' vectors.
#'
#' @details
#' The Kaplan-Meier estimate at time \code{t_eval} is defined as the
#' product-limit estimator evaluated at the largest observed event time less
#' than or equal to \code{t_eval}. If \code{t_eval} is smaller than the first
#' observed event time, \code{S(t) = 1} and the standard error is zero.
#'
#' The standard error is estimated by Greenwood's formula:
#'
#' SE[S(t)] = S(t) * sqrt(sum_\{t_i <= t, d_i > 0\} d_i / (n_i * (n_i - d_i)))
#'
#' where d_i is the number of events and n_i is the number at risk at time t_i.
#' The output field \code{std.err} follows the convention of
#' \code{\link[survival]{survfit}}, which reports SE[S(t)] / S(t) (i.e., the
#' standard error on the log scale) when \code{conf.type != "plain"}, and
#' SE[S(t)] when \code{conf.type = "plain"}. This function always returns
#' SE[S(t)] (the standard error on the survival scale).
#'
#' When \code{S(t_eval) = 0} (all subjects have experienced the event by
#' \code{t_eval}), the standard error is zero and the confidence interval
#' collapses to \code{[0, 0]}, consistent with \code{\link[survival]{survfit}}.
#'
#' When \code{presorted = TRUE} (default), \code{t_sorted} and \code{e_sorted}
#' are assumed to be sorted in ascending order of time. When
#' \code{presorted = FALSE}, the vectors are sorted internally before
#' computation.
#'
#' Three confidence interval types are supported via \code{conf.type}:
#'
#' \itemize{
#'   \item \code{"plain"}: Linear interval on the survival scale,
#'     S(t) +/- z * SE. The bounds are clipped to [0, 1].
#'   \item \code{"log"}: Interval on the log scale (default in
#'     \code{\link[survival]{survfit}}),
#'     S(t) * exp(+/- z * SE / S(t)).
#'   \item \code{"log-log"}: Interval on the complementary log-log scale,
#'     S(t)^exp(+/- z * SE / (S(t) * log(S(t)))).
#' }
#'
#' The returned object has class \code{"survfit_fast"} and is a named numeric
#' vector of length 4 with the evaluation time \code{t_eval}, the confidence
#' level \code{conf.level}, and the confidence interval type \code{conf.type}
#' stored as attributes. A \code{print()} method formats the result similarly
#' to \code{print(summary(survival::survfit(...)))}.
#'
#' @param t_sorted A numeric vector of event or censoring times. Must be
#'   sorted in ascending order when \code{presorted = TRUE}.
#' @param e_sorted An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{t_sorted}.
#' @param t_eval A single numeric value specifying the time point at which
#'   the survival probability is evaluated.
#' @param conf.level A single numeric value in (0, 1) specifying the
#'   confidence level. Defaults to 0.95.
#' @param conf.type A character string specifying the confidence interval type.
#'   Must be one of \code{"plain"}, \code{"log"}, or \code{"log-log"}.
#'   Defaults to \code{"log"}.
#' @param presorted A logical value. If \code{TRUE} (default), \code{t_sorted}
#'   and \code{e_sorted} are assumed to be sorted in ascending order of time.
#'   If \code{FALSE}, the vectors are sorted internally before computation.
#'
#' @return An object of class \code{"survfit_fast"}, which is a named numeric
#'   vector of length 4 with elements \code{surv}, \code{std.err},
#'   \code{lower}, and \code{upper}, representing the Kaplan-Meier survival
#'   estimate, the Greenwood standard error SE[S(t)], and the lower and upper
#'   confidence limits at \code{t_eval}. The evaluation time, confidence
#'   level, and confidence interval type are stored as attributes
#'   \code{t_eval}, \code{conf.level}, and \code{conf.type}.
#'   Returns a vector of \code{NA_real_} values (still with class
#'   \code{"survfit_fast"}) when \code{n} is zero.
#'
#' @examples
#' set.seed(42)
#' t_raw <- rexp(100, rate = 1 / 10)
#' e_raw <- rbinom(100, 1, 0.7)
#'
#' # presorted = TRUE (default): sort once outside, reuse inside a loop
#' ord <- order(t_raw)
#' t_s <- t_raw[ord]
#' e_s <- e_raw[ord]
#' survfit_fast(t_s, e_s, t_eval = 10, conf.type = "plain")
#' survfit_fast(t_s, e_s, t_eval = 10, conf.type = "log")
#' survfit_fast(t_s, e_s, t_eval = 10, conf.type = "log-log")
#'
#' # presorted = FALSE: sort internally, convenient for one-off calls
#' survfit_fast(t_raw, e_raw, t_eval = 10, presorted = FALSE)
#'
#' # Validation against survival::survfit
#' library(survival)
#' fit <- survfit(Surv(t_raw, e_raw) ~ 1, conf.type = "plain")
#' summary(fit, times = 10)
#'
#' @seealso
#' \code{\link[survival]{survfit}} for the standard Kaplan-Meier estimator.
#' \code{\link{print.survfit_fast}} for the print method.
#'
#' @references
#' Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete
#' observations. \emph{Journal of the American Statistical Association},
#' \emph{53}(282), 457-481.
#'
#' @importFrom stats qnorm
#' @export
survfit_fast <- function(t_sorted, e_sorted, t_eval,
                         conf.level = 0.95, conf.type = "log",
                         presorted = TRUE) {
  conf.type <- match.arg(conf.type, choices = c("plain", "log", "log-log"))
  na_out <- c(surv = NA_real_, std.err = NA_real_,
              lower = NA_real_, upper = NA_real_)
  wrap <- function(v) {
    structure(v,
              t_eval    = t_eval,
              conf.level  = conf.level,
              conf.type = conf.type,
              class     = "survfit_fast")
  }

  n <- length(t_sorted)
  if (n == 0L) return(wrap(na_out))

  # Sort internally when presorted = FALSE
  if (!presorted) {
    ord      <- order(t_sorted)
    t_sorted <- t_sorted[ord]
    e_sorted <- e_sorted[ord]
  }

  # C++ core: binary search + single scan -> c(surv, gw_sum)
  res    <- km_core(t_sorted, e_sorted, t_eval)
  surv   <- res[1L]
  gw_sum <- res[2L]

  if (is.na(surv)) return(wrap(na_out))
  if (surv == 0) {
    return(wrap(c(surv = 0, std.err = 0, lower = 0, upper = 0)))
  }

  std.err <- surv * sqrt(gw_sum)
  z       <- qnorm(1 - (1 - conf.level) / 2)

  ci <- if (conf.type == "plain") {
    c(lower = max(0, surv - z * std.err),
      upper = min(1, surv + z * std.err))
  } else if (conf.type == "log") {
    c(lower = surv * exp(-z * std.err / surv),
      upper = surv * exp( z * std.err / surv))
  } else {
    if (surv >= 1) {
      return(wrap(c(surv = surv, std.err = std.err,
                    lower = NA_real_, upper = NA_real_)))
    }
    log_s    <- log(surv)
    se_ll    <- std.err / (surv * abs(log_s))
    theta_lo <- log(-log_s) - z * se_ll
    theta_hi <- log(-log_s) + z * se_ll
    c(lower = exp(-exp(theta_hi)),
      upper = exp(-exp(theta_lo)))
  }

  wrap(c(surv = surv, std.err = std.err, ci))
}
