#' Print Method for wmst_fast Objects
#'
#' @description
#' Formats and prints a \code{wmst_fast} object in the same layout as the other
#' two-group summaries in the package. The header shows the window, the control
#' label and the inference settings. The body shows the per-group window mean
#' survival time with its confidence interval, followed, for a two-group object,
#' by the difference contrast (treatment minus control) with a confidence
#' interval, the test statistic, and the p-value.
#'
#' @param x An object of class \code{"wmst_fast"} returned by
#'   \code{\link{wmst_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' g <- rep(0:1, each = n / 2)
#' tt <- c(rexp(n / 2, 0.1), rexp(n / 2, 0.07))
#' cc <- rexp(n, 0.02)
#' time <- pmin(tt, cc)
#' event <- as.integer(tt <= cc)
#' print(wmst_fast(time, event, group = g, control = 0, tau1 = 2, tau2 = 12))
#'
#' @seealso \code{\link{wmst_fast}}
#'
#' @importFrom stats printCoefmat
#' @export
#' @method print wmst_fast
print.wmst_fast <- function(x, digits = max(1L, getOption("digits") - 3L),
                            ...) {
  two_group <- attr(x, "two.group")
  conf.level <- attr(x, "conf.level")
  if (is.null(conf.level)) conf.level <- 0.95
  ci_lab <- conf.level * 100
  tau1 <- attr(x, "tau1")
  tau2 <- attr(x, "tau2")
  side <- attr(x, "side")
  alt <- if (side == 2) "two.sided" else "one.sided"

  if (two_group) {
    control <- attr(x, "control")
    cat("Window mean survival time (two-group)\n\n")
    cat(sprintf("  window = [%g, %g],  control = %s\n",
                tau1, tau2, format(control)))
    cat(sprintf("  alternative = %s\n\n", alt))

    g_tab <- matrix(
      as.numeric(c(x["wmst.control"], x["wmst.treatment"],
                   x["se.control"], x["se.treatment"],
                   x["lower.control"], x["lower.treatment"],
                   x["upper.control"], x["upper.treatment"])),
      nrow = 2L,
      dimnames = list(
        c("control", "treatment"),
        c("wmst", "std.err",
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
    cat("Window mean survival time (single-group)\n\n")
    cat(sprintf("  window = [%g, %g]\n\n", tau1, tau2))

    g_tab <- matrix(
      as.numeric(c(x["wmst"], x["se"], x["lower"], x["upper"])),
      nrow = 1L,
      dimnames = list(
        "group",
        c("wmst", "std.err",
          sprintf("lower %g%%", ci_lab),
          sprintf("upper %g%%", ci_lab))
      )
    )
    print(round(g_tab, digits = digits))
  }

  invisible(x)
}
