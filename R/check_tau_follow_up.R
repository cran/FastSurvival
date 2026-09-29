#' Warn when a truncation time exceeds the observed follow-up
#'
#' Internal helper that warns when \code{tau} is larger than the largest
#' observed time (event or censoring) in a group. Beyond that time the
#' Kaplan-Meier curve of the group is not estimated and is carried forward
#' flat, which reference implementations such as survRM2 refuse. For two
#' groups the limit is the smaller of the two group maxima. Nothing is checked
#' when a group is empty.
#'
#' @param time A numeric vector of follow-up times without missing values.
#' @param j An integer vector of group indicators (0 for control, 1 for
#'   treatment), or \code{NULL} for a single group.
#' @param tau A single numeric value, the truncation time or milestone.
#' @param arg The argument name used in the warning.
#'
#' @return \code{NULL}, invisibly.
#'
#' @keywords internal
#' @noRd
check_tau_follow_up <- function(time, j, tau, arg = "tau") {
  if (is.null(j)) {
    if (length(time) == 0L) return(invisible(NULL))
    limit <- max(time)
    where <- "the largest observed time"
  } else {
    t0 <- time[j == 0L]
    t1 <- time[j == 1L]
    if (length(t0) == 0L || length(t1) == 0L) return(invisible(NULL))
    limit <- min(max(t0), max(t1))
    where <- "the largest observed time of one group"
  }
  if (tau > limit) {
    warning("'", arg, "' (", format(tau), ") exceeds ", where, " (",
            format(limit), "); the Kaplan-Meier curve is carried forward ",
            "flat beyond that time", call. = FALSE)
  }
  invisible(NULL)
}
