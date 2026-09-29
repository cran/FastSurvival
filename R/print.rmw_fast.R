#' Print Method for rmw_fast Objects
#'
#' @description
#' Formats and prints an \code{rmw_fast} object, showing the standardized
#' log-rank and modestly-weighted component Z-scores, their null correlation,
#' the survival-probability threshold \code{s_star}, the combined test
#' statistic, and the corresponding p-value.
#'
#' @param x An object of class \code{"rmw_fast"} returned by
#'   \code{\link{rmw_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' library(survival)
#' fit <- rmw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
#'                 control = 1, side = 1, s_star = 0.5)
#' print(fit)
#'
#' @seealso \code{\link{rmw_fast}}
#'
#' @export
#' @method print rmw_fast
print.rmw_fast <- function(x, digits = max(1L, getOption("digits") - 3L), ...) {

  side   <- attr(x, "side")
  z      <- attr(x, "z")
  corr   <- attr(x, "corr")
  s_star <- attr(x, "s_star")
  n      <- attr(x, "n")
  control <- attr(x, "control")

  cat("Robust modestly-weighted log-rank test (two-group)\n\n")
  if (is.null(control)) {
    cat(sprintf("  N = %d,  s_star = %s\n", n, format(s_star)))
  } else {
    cat(sprintf("  N = %d,  control = %s,  s_star = %s\n",
                n, format(control), format(s_star)))
  }

  stat_val <- as.numeric(x)[1L]
  p_val    <- as.numeric(x)[2L]

  if (is.na(stat_val)) {
    cat("  Test statistic not available (zero variance).\n")
    return(invisible(x))
  }

  ztab <- matrix(
    c(z[["logrank"]], z[["mwlrt"]]),
    nrow = 2L,
    dimnames = list(c("log-rank", "modestly-weighted"), "Z")
  )
  cat("\n")
  print(round(ztab, digits = digits))
  cat(sprintf("\n  Null correlation = %s\n", format(round(corr[1L, 2L], digits))))

  if (side == 2L) {
    cat(sprintf(" Max |Z| = %s,  two-sided p-value = %s %s\n",
                format(signif(stat_val, digits)),
                format.pval(p_val, digits = digits), signif_star(p_val)))
  } else {
    cat(sprintf(" min Z = %s,  one-sided p-value = %s %s\n",
                format(signif(stat_val, digits)),
                format.pval(p_val, digits = digits), signif_star(p_val)))
  }
  cat(signif_legend())

  invisible(x)
}
