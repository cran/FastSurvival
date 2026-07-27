#' Print Method for wkm_fast Objects
#'
#' @description
#' Formats and prints a \code{wkm_fast} object in the same layout as the other
#' two-group summaries in the package. The header shows the control label and
#' the inference settings. The body shows the weighted integrated survival
#' difference (treatment minus control) with a confidence interval, the test
#' statistic, and the p-value.
#'
#' @param x An object of class \code{"wkm_fast"} returned by
#'   \code{\link{wkm_fast}}.
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
#' print(wkm_fast(time, event, group = g, control = 0))
#'
#' @seealso \code{\link{wkm_fast}}
#'
#' @importFrom stats printCoefmat
#' @export
#' @method print wkm_fast
print.wkm_fast <- function(x, digits = max(1L, getOption("digits") - 3L),
                           ...) {
  conf.level <- attr(x, "conf.level")
  if (is.null(conf.level)) conf.level <- 0.95
  ci_lab <- conf.level * 100
  control <- attr(x, "control")
  weight <- attr(x, "weight")
  side <- attr(x, "side")
  alt <- if (side == 2) "two.sided" else "one.sided"

  cat("Weighted Kaplan-Meier test (Pepe-Fleming, two-group)\n\n")
  cat(sprintf("  control = %s\n", format(control)))
  cat(sprintf("  weight = %s,  alternative = %s\n\n", weight, alt))

  if (is.na(x["z"])) {
    cat("  Test statistic not available.\n")
    return(invisible(x))
  }

  con_mat <- matrix(
    as.numeric(c(x["wdiff"], x["lower"], x["upper"], x["z"], x["p"])),
    nrow = 1L,
    dimnames = list(
      "weighted difference (treatment - control)",
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
