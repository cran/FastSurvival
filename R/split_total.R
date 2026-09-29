#' Split a total sample size between two groups
#'
#' Internal helper that divides a total sample size between the control and
#' treatment groups in proportion to an allocation ratio. Each group receives
#' the integer part of its share, and the remaining subjects go to the groups
#' with the largest fractional parts, so the two sizes always add up to
#' \code{n}.
#'
#' @param n A single positive whole number, the total sample size.
#' @param alloc A numeric vector of length two with positive entries, the
#'   allocation ratio (control first, treatment second).
#'
#' @return An integer vector of length two, the control and treatment sizes.
#'
#' @keywords internal
#' @noRd
split_total <- function(n, alloc) {
  if (length(n) != 1L || !is.finite(n) || n < 1 ||
      abs(n - round(n)) > 1e-8) {
    stop("A scalar 'n' must be a positive whole number")
  }
  if (!is.numeric(alloc) || length(alloc) != 2L || any(!is.finite(alloc)) ||
      any(alloc <= 0)) {
    stop("'alloc' must be a numeric vector of two positive values")
  }
  accrual_cell_counts(as.integer(round(n)), alloc / sum(alloc))
}
