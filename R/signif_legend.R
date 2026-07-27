#' Significance code legend
#'
#' Returns the significance-code legend line used by \code{stats::printCoefmat},
#' so that the test summaries printed with \code{cat()} carry the same legend as
#' the contrast tables printed with \code{printCoefmat}.
#'
#' @return A length-one character string ending in a newline.
#' @noRd
signif_legend <- function() {
  "---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n"
}
