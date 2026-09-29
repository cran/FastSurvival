#' Print Method for milestone_fast Objects
#'
#' @description
#' Formats and prints a \code{milestone_fast} object in the same layout as the
#' other two-group estimation summaries in the package. The header shows the
#' milestone timepoint, the control label, and the inference settings. The body
#' shows the per-group milestone survival with its confidence interval, followed
#' by the difference contrast (treatment minus control) with a confidence
#' interval, the test statistic, and the p-value. The p-value follows the
#' alternative recorded in the object.
#'
#' @param x An object of class \code{"milestone_fast"} returned by
#'   \code{\link{milestone_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(1)
#' time <- c(rexp(50, 0.1), rexp(50, 0.07))
#' status <- rep(1, 100)
#' group <- rep(c(0, 1), each = 50)
#' print(milestone_fast(time, status, group, control = 0, tau = 10, method = "loglog"))
#'
#' @seealso \code{\link{milestone_fast}}
#'
#' @importFrom stats printCoefmat
#' @export
#' @method print milestone_fast
print.milestone_fast <- function(x, digits = max(1L, getOption("digits") - 3L),
                                 ...) {
  tau        <- x[["tau"]]
  conf.level <- x[["conf.level"]]
  if (is.null(conf.level)) conf.level <- 0.95
  ci_lab     <- conf.level * 100
  control    <- x[["group.labels"]][["control"]]

  cat("Milestone survival (two-group)\n\n")
  cat(sprintf("  tau = %g,  control = %s\n", tau, format(control)))
  alt <- x[["side"]]
  if (is.numeric(alt)) alt <- if (alt == 1L) "one.sided" else "two.sided"
  cat(sprintf("  method = %s,  alternative = %s\n\n",
              x[["method"]], alt))

  if (is.na(x[["surv"]][["control"]]) || is.na(x[["surv"]][["treatment"]])) {
    cat("  Estimate not available.\n")
    return(invisible(x))
  }

  # Per-group milestone survival.
  g_tab <- matrix(
    c(x[["surv"]], x[["std.err"]], x[["surv.lower"]], x[["surv.upper"]]),
    nrow = 2L,
    dimnames = list(
      c("control", "treatment"),
      c("survival", "std.err",
        sprintf("lower %g%%", ci_lab),
        sprintf("upper %g%%", ci_lab))
    )
  )
  print(round(g_tab, digits = digits))
  cat("\n")

  # Difference contrast (treatment minus control).
  con_mat <- matrix(
    c(x[["diff"]], x[["diff.lower"]], x[["diff.upper"]],
      x[["statistic"]], x[["p.value"]]),
    nrow = 1L,
    dimnames = list(
      "difference (treatment - control)",
      c("Est.",
        sprintf("lower %g%%", ci_lab),
        sprintf("upper %g%%", ci_lab),
        "z", "Pr(>|z|)")
    )
  )
  printCoefmat(con_mat, digits = digits, has.Pvalue = TRUE, P.values = TRUE,
               cs.ind = 1:3, tst.ind = 4L)

  invisible(x)
}
