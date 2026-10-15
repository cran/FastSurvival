#' Numbers at risk of a kmcurve_fast object
#'
#' Internal helper of \code{plot.kmcurve_fast()} that counts, for each group,
#' the subjects whose observed time is at or after each of the given times.
#'
#' @param x A \code{kmcurve_fast} object.
#' @param times A numeric vector of times.
#'
#' @return A numeric matrix with one row per group (\code{control}, and
#'   \code{treat} in the two-group case) and one column per time.
#'
#' @keywords internal
#' @noRd
km_risk_counts <- function(x, times) {
  dat   <- x$data
  count <- function(tt) {
    vapply(times, function(u) as.numeric(sum(tt >= u)), numeric(1L))
  }
  if (isTRUE(x$two_group)) {
    is_ctrl <- as.character(dat$group) == as.character(dat$control)
    rbind(control = count(dat$time[is_ctrl]),
          treat   = count(dat$time[!is_ctrl]))
  } else {
    matrix(count(dat$time), nrow = 1L, dimnames = list("control", NULL))
  }
}
