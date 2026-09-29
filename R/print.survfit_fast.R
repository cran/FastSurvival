#' Print Method for survfit_fast Objects
#'
#' @description
#' Formats and prints a \code{survfit_fast} object similarly to
#' \code{print(summary(survival::survfit(...), times = t_eval))}, showing the
#' Kaplan-Meier survival estimate, the Greenwood standard error on the
#' survival scale, and the confidence interval at the requested evaluation
#' time.
#'
#' @param x An object of class \code{"survfit_fast"} returned by
#'   \code{\link{survfit_fast}}.
#' @param digits Number of significant digits to display. Defaults to the
#'   global option \code{getOption("digits")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(42)
#' t_raw <- rexp(100, rate = 1 / 10)
#' e_raw <- rbinom(100, 1, 0.7)
#' ord   <- order(t_raw)
#' fit   <- survfit_fast(t_raw[ord], e_raw[ord], t_eval = 10)
#' print(fit)
#'
#' @seealso \code{\link{survfit_fast}}
#'
#' @export
#' @method print survfit_fast
print.survfit_fast <- function(x, digits = max(1L, getOption("digits") - 3L), ...) {

  t_eval    <- attr(x, "t_eval")
  conf.level  <- attr(x, "conf.level")
  conf.type <- attr(x, "conf.type")
  if (is.null(conf.level))  conf.level  <- 0.95
  if (is.null(conf.type)) conf.type <- "log"
  ci_lab <- conf.level * 100

  cat("Kaplan-Meier survival estimate (single time point)\n\n")

  if (anyNA(x)) {
    cat("  Estimate not available.\n")
    return(invisible(x))
  }

  surv_v <- x[["surv"]]
  se_v   <- x[["std.err"]]
  lo_v   <- x[["lower"]]
  hi_v   <- x[["upper"]]

  tab <- matrix(
    c(surv_v, se_v, lo_v, hi_v),
    nrow = 1L,
    dimnames = list(
      sprintf("t = %g", t_eval),
      c("survival", "std.err",
        sprintf("lower %g%%", ci_lab),
        sprintf("upper %g%%", ci_lab))
    )
  )

  print(round(tab, digits = digits))
  cat(sprintf("\n Confidence interval type: %s\n", conf.type))

  invisible(x)
}
