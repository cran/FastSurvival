#' Significance star for a single p-value
#'
#' Returns the significance code for a single p-value using the same cutpoints
#' and symbols as \code{stats::printCoefmat}, so that the test summaries printed
#' with \code{cat()} carry the same stars as the contrast tables printed with
#' \code{printCoefmat}. Returns an empty string for a missing or non-finite
#' value.
#'
#' @param p A single numeric p-value.
#' @return A length-one character string, one of \code{"***"}, \code{"**"},
#'   \code{"*"}, \code{"."}, \code{" "}, or \code{""}.
#' @noRd
signif_star <- function(p) {
  if (length(p) != 1L || is.na(p) || !is.finite(p)) {
    return("")
  }
  as.character(stats::symnum(p, corr = FALSE, na = FALSE,
                             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " ")))
}
