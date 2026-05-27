#' Print Method for coxph_fast Objects
#'
#' @description
#' Formats and prints a \code{coxph_fast} object similarly to
#' \code{summary(survival::coxph(...))}, showing the point estimate of the
#' log hazard ratio, the hazard ratio, the standard error on the log scale,
#' the Wald z-statistic, the corresponding two-sided p-value, and the Wald
#' confidence interval for the hazard ratio.
#'
#' @param x An object of class \code{"coxph_fast"} returned by
#'   \code{\link{coxph_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments passed to \code{\link[base]{format}}
#'   (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' library(survival)
#' fit <- coxph_fast(ovarian$futime, ovarian$fustat, ovarian$rx, control = 1)
#' print(fit)
#'
#' @seealso \code{\link{coxph_fast}}
#'
#' @importFrom stats pnorm printCoefmat
#' @export
#' @method print coxph_fast
print.coxph_fast <- function(x, digits = max(1L, getOption("digits") - 3L), ...) {

  # Detect whether the estimate is available
  if (anyNA(x)) {
    cat("Pike-Halley Estimator (closed-form hazard ratio)\n\n")
    cat("  Estimate not available (insufficient data).\n")
    return(invisible(x))
  }

  coef_v    <- x[["coef"]]
  hr_v      <- x[["exp(coef)"]]
  se_v      <- x[["se(coef)"]]
  conf.int  <- attr(x, "conf.int")
  if (is.null(conf.int)) conf.int <- 0.95
  ci_lab    <- conf.int * 100
  lo_v      <- x[[sprintf("lower .%g", ci_lab)]]
  hi_v      <- x[[sprintf("upper .%g", ci_lab)]]

  # Wald z and two-sided p-value
  z_v <- coef_v / se_v
  p_v <- 2 * pnorm(-abs(z_v))

  # Coefficient table: matches the column layout of summary(coxph(...))$coefficients
  coef_mat <- matrix(
    c(coef_v, hr_v, se_v, z_v, p_v),
    nrow = 1L,
    dimnames = list(
      "group",
      c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
    )
  )

  # Confidence interval table: matches summary(coxph(...))$conf.int
  ci_mat <- matrix(
    c(hr_v, 1 / hr_v, lo_v, hi_v),
    nrow = 1L,
    dimnames = list(
      "group",
      c("exp(coef)", "exp(-coef)",
        sprintf("lower .%g", ci_lab),
        sprintf("upper .%g", ci_lab))
    )
  )

  cat("Pike-Halley Estimator (closed-form hazard ratio)\n\n")
  cat("Coefficients:\n")
  printCoefmat(coef_mat, digits = digits, has.Pvalue = TRUE, P.values = TRUE)
  cat("\n")
  cat(sprintf("Hazard ratio and %g%% Wald confidence interval:\n", ci_lab))
  print(round(ci_mat, digits = digits))

  invisible(x)
}
