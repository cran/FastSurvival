#' Median of a Kaplan-Meier step function
#'
#' Internal helper that reads the median from a Kaplan-Meier curve with the
#' convention of \code{survival::survfit}: the first event time at which the
#' survival estimate is 0.5 or below, and, when the estimate equals 0.5 on a
#' flat stretch, the midpoint between that event time and the next event time
#' (or \code{tmax} when there is no later event). The comparison with 0.5 uses
#' the tolerance \code{sqrt(.Machine$double.eps)}. When that stretch ends at an
#' infinite time, the midpoint is not defined and \code{NA} is returned.
#'
#' @param te A numeric vector of distinct event times in ascending order.
#' @param surv A numeric vector of Kaplan-Meier estimates at \code{te}.
#' @param tmax The largest observed time, used to close a final flat stretch.
#'
#' @return A single numeric value, \code{NA} when the curve does not reach
#'   0.5 or stays at 0.5 up to an infinite time.
#'
#' @keywords internal
#' @noRd
km_step_median <- function(te, surv, tmax) {
  tol <- sqrt(.Machine$double.eps)
  if (length(te) == 0L) return(NA_real_)
  i <- which(surv <= 0.5 + tol)
  if (length(i) == 0L) return(NA_real_)
  i <- i[1L]
  if (abs(surv[i] - 0.5) < tol) {
    t_next <- if (i < length(te)) te[i + 1L] else tmax
    if (!is.finite(t_next)) return(NA_real_)
    return((te[i] + t_next) / 2)
  }
  te[i]
}
