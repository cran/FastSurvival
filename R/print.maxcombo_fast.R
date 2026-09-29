#' Print Method for maxcombo_fast Objects
#'
#' @description
#' Formats and prints a \code{maxcombo_fast} object. It shows the total sample
#' size, a table of the Fleming-Harrington component Z-scores, and the
#' max-combo statistic with its p-value. For a one-sided test the statistic is
#' the most negative component, so a negative value favors the treatment
#' group, matching the sign convention of \code{\link{survdiff_fast}}. For a
#' two-sided test the statistic is the largest component in absolute value.
#'
#' @param x An object of class \code{"maxcombo_fast"} returned by
#'   \code{\link{maxcombo_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' library(survival)
#' fit <- maxcombo_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
#'                      control = 1, side = 1)
#' print(fit)
#'
#' @seealso \code{\link{maxcombo_fast}}
#'
#' @export
#' @method print maxcombo_fast
print.maxcombo_fast <- function(x, digits = max(1L, getOption("digits") - 3L),
                                ...) {
  side  <- attr(x, "side")
  z     <- attr(x, "z")
  rho   <- attr(x, "rho")
  gamma <- attr(x, "gamma")
  n     <- attr(x, "n")
  control <- attr(x, "control")

  cat("Max-combo weighted log-rank test (two-group)\n\n")
  if (is.null(control)) {
    cat(sprintf("  N = %d\n", n))
  } else {
    cat(sprintf("  N = %d,  control = %s\n", n, format(control)))
  }

  stat_val <- as.numeric(x["statistic"])
  p_val    <- as.numeric(x["p.value"])

  if (is.na(stat_val)) {
    cat("  Test statistic not available (a component variance is zero).\n")
    return(invisible(x))
  }

  # Table of component Z-scores, one row per Fleming-Harrington weight.
  tab <- matrix(
    z,
    ncol = 1L,
    dimnames = list(sprintf("FH(%g,%g)", rho, gamma), "Z")
  )
  cat("\n")
  print(round(tab, digits = digits))
  cat("\n")

  if (side == 1L) {
    cat(sprintf(" Max-combo statistic = %s (one-sided),  p-value = %s %s\n",
                format(signif(stat_val, digits)),
                format.pval(p_val, digits = digits), signif_star(p_val)))
  } else {
    cat(sprintf(" Max-combo statistic = %s (two-sided),  p-value = %s %s\n",
                format(signif(stat_val, digits)),
                format.pval(p_val, digits = digits), signif_star(p_val)))
  }
  cat(signif_legend())

  invisible(x)
}
