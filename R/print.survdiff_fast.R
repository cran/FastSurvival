#' Print Method for survdiff_fast Objects
#'
#' @description
#' Formats and prints a \code{survdiff_fast} object similarly to
#' \code{print(survival::survdiff(...))}, showing the observed and expected
#' event counts for the control and treatment groups, the per-group
#' contributions \code{(O-E)^2 / E} and \code{(O-E)^2 / V}, the test
#' statistic, and the corresponding p-value.
#'
#' @param x An object of class \code{"survdiff_fast"} returned by
#'   \code{\link{survdiff_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' library(survival)
#' fit <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
#'                      control = 1, side = 2)
#' print(fit)
#'
#' @seealso \code{\link{survdiff_fast}}
#'
#' @importFrom stats pchisq pnorm
#' @export
#' @method print survdiff_fast
print.survdiff_fast <- function(x, digits = max(1L, getOption("digits") - 3L), ...) {

  side <- attr(x, "side")
  O0   <- attr(x, "O0")
  E0   <- attr(x, "E0")
  O1   <- attr(x, "O1")
  E1   <- attr(x, "E1")
  V1   <- attr(x, "V1")
  n    <- attr(x, "n")

  cat("Log-rank test (two-group)\n\n")
  cat(sprintf("  N = %d\n", n))

  if (is.na(x) || !is.finite(V1) || V1 == 0) {
    cat("  Test statistic not available (V1 = 0).\n")
    return(invisible(x))
  }

  # Two-row table: control (row 0) and treatment (row 1).
  # By the log-rank identity (O0 - E0) = -(O1 - E1), so (O0 - E0)^2 = (O1 - E1)^2
  # and the per-group (O - E)^2 / V contribution is identical for both rows.
  OE_sq <- (O1 - E1) ^ 2
  tab <- matrix(
    c(O0, O1,
      E0, E1,
      OE_sq / E0, OE_sq / E1,
      OE_sq / V1, OE_sq / V1),
    nrow = 2L,
    dimnames = list(
      c("control", "treatment"),
      c("Observed", "Expected", "(O-E)^2/E", "(O-E)^2/V")
    )
  )
  cat("\n")
  print(round(tab, digits = digits))
  cat("\n")

  stat_val <- as.numeric(x)
  if (side == 2L) {
    p_val <- pchisq(stat_val, df = 1L, lower.tail = FALSE)
    cat(sprintf(" Chi-square = %s on 1 df,  p-value = %s\n",
                format(signif(stat_val, digits)),
                format.pval(p_val, digits = digits)))
  } else {
    p_val <- 2 * pnorm(-abs(stat_val))
    cat(sprintf(" Z = %s,  two-sided p-value = %s\n",
                format(signif(stat_val, digits)),
                format.pval(p_val, digits = digits)))
  }

  invisible(x)
}
