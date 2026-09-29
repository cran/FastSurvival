#' Print Method for ahsw_fast Objects
#'
#' @description
#' Formats and prints an \code{ahsw_fast} object. The header shows the
#' truncation time and the control label. The body shows the per-group average
#' hazard with survival weight, followed by the between-group contrasts: the
#' ratio of average hazards (treatment over control) and the difference of
#' average hazards (treatment minus control), each with a confidence interval
#' and a two-sided p-value.
#'
#' @param x An object of class \code{"ahsw_fast"} returned by
#'   \code{\link{ahsw_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' library(survival)
#' fit <- ahsw_fast(ovarian$futime, ovarian$fustat, ovarian$rx,
#'                  control = 1, tau = 600)
#' print(fit)
#'
#' @seealso \code{\link{ahsw_fast}}
#'
#' @importFrom stats printCoefmat
#' @export
#' @method print ahsw_fast
print.ahsw_fast <- function(x, digits = max(1L, getOption("digits") - 3L),
                            ...) {
  tau      <- attr(x, "tau")
  conf.level <- attr(x, "conf.level")
  if (is.null(conf.level)) conf.level <- 0.95
  ci_lab   <- conf.level * 100
  control  <- attr(x, "control")
  side     <- attr(x, "side")
  if (is.null(side)) side <- 2L

  cat("Average hazard with survival weight (two-group)\n\n")
  cat(sprintf("  tau = %g,  control = %s\n", tau, format(control)))
  cat(sprintf("  alternative = %s\n\n",
              if (side == 1L) "one.sided" else "two.sided"))

  if (is.na(x[["ah.ctrl"]]) || is.na(x[["ah.trt"]])) {
    cat("  Estimate not available (no events up to tau in a group or",
        "non-finite variance).\n")
    return(invisible(x))
  }

  # Per-group average hazard.
  g_tab <- matrix(
    c(x[["ah.ctrl"]], x[["ah.trt"]]),
    nrow = 2L,
    dimnames = list(c("control", "treatment"), "AH")
  )
  print(round(g_tab, digits = digits))
  cat("\n")

  # Between-group contrasts: ratio (RAH) and difference (DAH).
  con_mat <- matrix(
    c(x[["rah"]], x[["rah.lower"]], x[["rah.upper"]], x[["p.rah"]],
      x[["dah"]], x[["dah.lower"]], x[["dah.upper"]], x[["p.dah"]]),
    nrow = 2L, byrow = TRUE,
    dimnames = list(
      c("ratio (treatment / control)", "difference (treatment - control)"),
      c("Est.",
        sprintf("lower %g%%", ci_lab),
        sprintf("upper %g%%", ci_lab),
        "p.value")
    )
  )
  printCoefmat(con_mat, digits = digits, P.values = TRUE, has.Pvalue = TRUE,
               cs.ind = 1:3, tst.ind = integer(0), zap.ind = integer(0))

  invisible(x)
}
