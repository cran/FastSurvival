#' Validate time and event vectors
#'
#' Internal helper that checks the follow-up times and event indicators passed
#' to the analysis functions. The times must be numeric, non-negative, and
#' without missing values, and the event indicator must have the same length and be coded as 0
#' (censored) or 1 (event) without missing values.
#'
#' @param time A numeric vector of follow-up times.
#' @param event A vector of event indicators.
#'
#' @return \code{NULL}, invisibly. An error is raised when a check fails.
#'
#' @keywords internal
#' @noRd
check_time_event <- function(time, event) {
  if (!is.numeric(time)) {
    stop("'time' must be numeric")
  }
  if (length(event) != length(time)) {
    stop("'time' and 'event' must have the same length")
  }
  if (anyNA(time)) {
    stop("'time' must not contain missing values")
  }
  if (any(time < 0)) {
    stop("'time' must be non-negative")
  }
  if (anyNA(event) || !all(event == 0 | event == 1)) {
    stop("'event' must be coded as 0 (censored) or 1 (event)")
  }
  invisible(NULL)
}
