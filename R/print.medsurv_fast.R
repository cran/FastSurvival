#' Print Method for medsurv_fast Objects
#'
#' @description
#' Formats and prints a \code{medsurv_fast} object in the same layout as the
#' other two-group estimation summaries in the package. The header shows the
#' control label and the inference settings. The body shows the per-group
#' median survival with its confidence interval, followed, for a two-group
#' object, by the difference contrast (treatment minus control) with a
#' confidence interval, the test statistic, and the p-value.
#'
#' @param x An object of class \code{"medsurv_fast"} returned by
#'   \code{\link{medsurv_fast}}.
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
#' print(medsurv_fast(time, status, group, control = 0))
#'
#' @seealso \code{\link{medsurv_fast}}
#'
#' @importFrom stats printCoefmat
#' @export
#' @method print medsurv_fast
print.medsurv_fast <- function(x, digits = max(1L, getOption("digits") - 3L),
                               ...) {
  two_group  <- attr(x, "two.group")
  conf.level <- attr(x, "conf.level")
  if (is.null(conf.level)) conf.level <- 0.95
  ci_lab <- conf.level * 100
  method <- attr(x, "method")
  side   <- attr(x, "side")
  alt    <- if (side == 2) "two.sided" else "one.sided"

  if (two_group) {
    control <- attr(x, "control")
    cat("Median survival time (two-group)\n\n")
    cat(sprintf("  control = %s\n", format(control)))
    cat(sprintf("  method = %s,  alternative = %s\n\n", method, alt))

    if (is.na(x["median.control"]) || is.na(x["median.treatment"])) {
      cat("  Estimate not available.\n")
      return(invisible(x))
    }

    g_tab <- matrix(
      as.numeric(c(x["median.control"], x["median.treatment"],
                   x["se.control"], x["se.treatment"],
                   x["lower.control"], x["lower.treatment"],
                   x["upper.control"], x["upper.treatment"])),
      nrow = 2L,
      dimnames = list(
        c("control", "treatment"),
        c("median", "std.err",
          sprintf("lower %g%%", ci_lab),
          sprintf("upper %g%%", ci_lab))
      )
    )
    print(round(g_tab, digits = digits))
    cat("\n")

    con_mat <- matrix(
      as.numeric(c(x["diff"], x["lower.diff"], x["upper.diff"],
                   x["z"], x["p"])),
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
  } else {
    cat("Median survival time (single-group)\n\n")
    cat(sprintf("  method = %s\n\n", method))

    if (is.na(x["median"])) {
      cat("  Estimate not available.\n")
      return(invisible(x))
    }

    g_tab <- matrix(
      as.numeric(c(x["median"], x["se"], x["lower"], x["upper"])),
      nrow = 1L,
      dimnames = list(
        "group",
        c("median", "std.err",
          sprintf("lower %g%%", ci_lab),
          sprintf("upper %g%%", ci_lab))
      )
    )
    print(round(g_tab, digits = digits))
  }

  invisible(x)
}
