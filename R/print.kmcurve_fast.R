#' Print Method for Kaplan-Meier Curves
#'
#' @description
#' Prints a compact per-group summary of a \code{kmcurve_fast} object: the group
#' label and role, the sample size, the number of events, and the median
#' survival read from the Kaplan-Meier curve. In the two-group case the constant
#' Cox proportional-hazards estimate of the hazard ratio is reported below the
#' table.
#'
#' @param x A \code{kmcurve_fast} object from \code{\link{kmcurve_fast}}.
#' @param ... Further arguments, currently ignored.
#'
#' @return The input \code{x}, invisibly.
#'
#' @seealso \code{\link{kmcurve_fast}}, \code{\link{plot.kmcurve_fast}}
#'
#' @export
print.kmcurve_fast <- function(x, ...) {
  two   <- isTRUE(x$two_group)
  km_c  <- x$km$control
  med_c <- curve_median(c(0, km_c$te), c(1, km_c$surv))

  if (two) {
    km_t  <- x$km$treat
    med_t <- curve_median(c(0, km_t$te), c(1, km_t$surv))
    tab <- data.frame(
      Group  = c(as.character(x$labels[["control"]]),
                 as.character(x$labels[["treat"]])),
      Role   = c("control", "treatment"),
      N      = c(km_c$n, km_t$n),
      Events = c(km_c$n_event, km_t$n_event),
      Median = round(c(med_c, med_t), 2),
      stringsAsFactors = FALSE
    )
  } else {
    tab <- data.frame(
      Group  = as.character(x$labels[["control"]]),
      Role   = "single",
      N      = km_c$n,
      Events = km_c$n_event,
      Median = round(med_c, 2),
      stringsAsFactors = FALSE
    )
  }

  cat(sprintf("A kmcurve_fast object (%s)\n\n",
              if (two) "two groups" else "one group"))
  print(tab, row.names = FALSE)

  if (two) {
    cx <- tryCatch(
      coxph_fast(x$data$time, x$data$event, x$data$group,
                 control = x$data$control),
      error = function(e) NULL)
    if (!is.null(cx) && is.finite(cx[["exp(coef)"]])) {
      cat(sprintf("\nCox PH hazard ratio (treatment vs control): %.3f\n",
                  unname(cx[["exp(coef)"]])))
    }
  }
  invisible(x)
}
