#' Abscissae of a Kaplan-Meier step curve
#'
#' Internal helper that returns time 0, the event times, and the largest
#' observed time when it lies beyond the last event time, so that a censored
#' tail of a Kaplan-Meier curve is drawn as a flat segment.
#'
#' @param te A numeric vector of distinct event times in ascending order.
#' @param tmax The largest observed time.
#'
#' @return A numeric vector of step abscissae.
#'
#' @keywords internal
#' @noRd
step_extend <- function(te, tmax) {
  xs <- c(0, te)
  if (length(tmax) == 1L && is.finite(tmax) && tmax > xs[length(xs)]) {
    xs <- c(xs, tmax)
  }
  xs
}
