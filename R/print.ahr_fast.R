#' Print an ahr_fast object
#'
#' @param x an object of class \code{"ahr_fast"}
#' @param digits number of significant digits to print
#' @param ... further arguments (currently ignored)
#' @return \code{x}, invisibly
#' @method print ahr_fast
#' @importFrom stats printCoefmat
#' @export
print.ahr_fast <- function(x, digits = max(3L, getOption("digits") - 3L),
                           ...) {
  ctrl <- as.character(x$groups[1])
  trt  <- as.character(x$groups[2])
  ci_pct <- round(100 * x$conf.level, 1)
  side <- if (is.null(x$side)) 2L else x$side
  alt  <- if (side == 1L) "one.sided" else "two.sided"

  cat("Kalbfleisch-Prentice average hazard ratio (two-group)\n\n")
  cat(sprintf("  tau = %s,  control = %s\n", format(x$tau, digits = digits),
              ctrl))
  cat(sprintf("  alternative = %s\n\n", alt))

  if (is.na(x$ahr)) {
    cat("  Estimate not available (insufficient data).\n")
    return(invisible(x))
  }

  # Per-group hazard shares, rows control and treatment.
  theta_ord <- c(x$theta[ctrl], x$theta[trt])
  n_ord     <- c(x$n[ctrl], x$n[trt])
  gmat <- cbind(theta = formatC(theta_ord, format = "f", digits = digits),
                n = formatC(n_ord, format = "d"))
  rownames(gmat) <- c("control", "treatment")
  print(gmat, quote = FALSE, right = TRUE)
  cat("\n")

  # Average hazard ratio contrast, treatment relative to control. The reported
  # test is the primary one on the theta (group-share) scale.
  cmat <- matrix(
    c(x$ahr, x$lower, x$upper, x$z, x$p.value),
    nrow = 1L,
    dimnames = list(
      "average hazard ratio (treatment / control)",
      c("Est.", sprintf("lower %g%%", ci_pct),
        sprintf("upper %g%%", ci_pct), "z", "Pr(>|z|)")
    )
  )
  printCoefmat(cmat, digits = digits, has.Pvalue = TRUE,
               P.values = TRUE, cs.ind = 1:3, tst.ind = 4L)
  cat(sprintf("  (log scale: z = %s, p = %s)\n",
              format(x$z.loghr, digits = digits),
              format.pval(x$p.value.loghr, digits = digits)))

  invisible(x)
}
