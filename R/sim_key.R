#' Character keys for simulation identifiers
#'
#' Internal helper that converts simulation identifiers to character keys used
#' as row names of per-simulation matrices (\code{\link{cutoff_fast}}) and for
#' matching them to \code{data$sim} (\code{\link{analysis_fast}},
#' \code{\link{switch_fast}}). Whole numbers are formatted without scientific
#' notation, so integer and double identifiers give the same keys.
#'
#' @param x A vector of simulation identifiers.
#'
#' @return A character vector of the same length as \code{x}.
#'
#' @keywords internal
#' @noRd
sim_key <- function(x) {
  if (is.numeric(x) && all(is.finite(x)) && all(x == round(x))) {
    format(round(x), scientific = FALSE, trim = TRUE)
  } else {
    as.character(x)
  }
}
