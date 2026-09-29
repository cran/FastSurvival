#' Print Method for rmst_fast Objects
#'
#' @description
#' Formats and prints an \code{rmst_fast} object. In single-group mode it shows
#' the restricted mean survival time, its Greenwood standard error, and the
#' Wald confidence interval at the requested horizon. In two-group mode it
#' shows the per-group restricted mean survival times together with the
#' difference (treatment minus control) and ratio (treatment over control)
#' contrasts, each with a Wald z-statistic and two-sided p-value.
#'
#' @param x An object of class \code{"rmst_fast"} returned by
#'   \code{\link{rmst_fast}}.
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
#'
#' # Single-group
#' print(rmst_fast(t_raw, e_raw, tau = 10))
#'
#' # Two-group comparison
#' set.seed(7)
#' n <- 200
#' time  <- c(rexp(n, 0.10), rexp(n, 0.07))
#' event <- rbinom(2 * n, 1, 0.8)
#' group <- rep(0:1, each = n)
#' print(rmst_fast(time, event, group = group, control = 0, tau = 10))
#'
#' @seealso \code{\link{rmst_fast}}
#'
#' @importFrom stats printCoefmat
#' @export
#' @method print rmst_fast
print.rmst_fast <- function(x, digits = max(1L, getOption("digits") - 3L), ...) {

  tau      <- attr(x, "tau")
  conf.level <- attr(x, "conf.level")
  if (is.null(conf.level)) conf.level <- 0.95
  ci_lab   <- conf.level * 100
  side     <- attr(x, "side")
  if (is.null(side)) side <- 2L

  two_group <- "diff" %in% names(x)

  # ---- Single-group mode -------------------------------------------------
  if (!two_group) {
    cat("Restricted mean survival time (single-group)\n\n")
    cat(sprintf("  tau = %g\n\n", tau))
    if (is.na(x[["rmst"]])) {
      cat("  Estimate not available.\n")
      return(invisible(x))
    }
    tab <- matrix(
      c(x[["rmst"]], x[["std.err"]], x[["lower"]], x[["upper"]]),
      nrow = 1L,
      dimnames = list(
        sprintf("tau = %g", tau),
        c("RMST", "std.err",
          sprintf("lower %g%%", ci_lab),
          sprintf("upper %g%%", ci_lab))
      )
    )
    print(round(tab, digits = digits))
    return(invisible(x))
  }

  # ---- Two-group mode ----------------------------------------------------
  control <- attr(x, "control")
  cat("Restricted mean survival time (two-group)\n\n")
  cat(sprintf("  tau = %g,  control = %s\n", tau, format(control)))
  cat(sprintf("  alternative = %s\n\n",
              if (side == 1L) "one.sided" else "two.sided"))

  g_tab <- matrix(
    c(x[["rmst.ctrl"]], x[["rmst.trt"]]),
    nrow = 2L,
    dimnames = list(c("control", "treatment"), "RMST")
  )
  print(round(g_tab, digits = digits))
  cat("\n")

  con_mat <- matrix(
    c(x[["diff"]],  x[["diff.lower"]],  x[["diff.upper"]],  x[["z.diff"]],  x[["p.diff"]],
      x[["ratio"]], x[["ratio.lower"]], x[["ratio.upper"]], x[["z.ratio"]], x[["p.ratio"]]),
    nrow = 2L, byrow = TRUE,
    dimnames = list(
      c("difference (treatment - control)", "ratio (treatment / control)"),
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
