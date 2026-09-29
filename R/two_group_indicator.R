#' Validate a two-group factor and build the treatment indicator
#'
#' Internal helper that checks that \code{group} has exactly two distinct
#' values without missing values and that \code{control} is one of them, and
#' returns the treatment indicator (1 for treatment, 0 for control).
#'
#' @param group A vector identifying the two groups. Factors are compared by
#'   their labels.
#' @param control The value of \code{group} that denotes the control group.
#'
#' @return An integer vector of the same length as \code{group}, 0 for the
#'   control group and 1 for the treatment group.
#'
#' @keywords internal
#' @noRd
two_group_indicator <- function(group, control) {
  if (missing(control) || is.null(control)) {
    stop("'control' must be supplied")
  }
  if (is.factor(group)) group <- as.character(group)
  if (anyNA(group)) {
    stop("'group' must not contain missing values")
  }
  lev <- unique(group)
  if (length(lev) != 2L) {
    stop("'group' must have exactly two distinct values")
  }
  if (length(control) != 1L || is.na(control) || !(control %in% lev)) {
    stop("'control' must be one of the two values in 'group'")
  }
  as.integer(group != control)
}
