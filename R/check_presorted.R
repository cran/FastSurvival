#' Check the order of presorted input
#'
#' Internal helper for the \code{presorted = TRUE} option of the analysis
#' functions. Without \code{strata}, \code{time} must be in ascending order.
#' With \code{strata}, the rows of each stratum must be contiguous and
#' \code{time} must be in ascending order within each stratum. The check is a
#' single pass over the data.
#'
#' @param time A numeric vector of follow-up times.
#' @param strata An optional vector of stratum labels, one per element of
#'   \code{time}.
#'
#' @return \code{NULL}, invisibly. An error is raised when the order does not
#'   hold.
#'
#' @keywords internal
#' @noRd
check_presorted <- function(time, strata = NULL) {
  n <- length(time)
  if (n < 2L) return(invisible(NULL))
  if (is.null(strata)) {
    if (is.unsorted(time)) {
      stop("with presorted = TRUE, 'time' must be sorted in ascending ",
           "order; use presorted = FALSE", call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (anyDuplicated(rle(as.vector(strata))$values) > 0L) {
    stop("with presorted = TRUE, the rows of each stratum must be ",
         "contiguous; use presorted = FALSE", call. = FALSE)
  }
  same <- strata[-1L] == strata[-n]
  if (any(same & time[-1L] < time[-n])) {
    stop("with presorted = TRUE, 'time' must be sorted in ascending order ",
         "within each stratum; use presorted = FALSE", call. = FALSE)
  }
  invisible(NULL)
}
