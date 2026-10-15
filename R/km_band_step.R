#' Shaded confidence band of a Kaplan-Meier curve
#'
#' Internal helper of \code{plot.kmcurve_fast()} that draws the pointwise
#' confidence limits of a Kaplan-Meier curve as a shaded band between the two
#' step functions, from time 0 to the largest observed time.
#'
#' @param km The Kaplan-Meier summary of one group of a \code{kmcurve_fast}
#'   object.
#' @param conf.level Confidence level of the limits.
#' @param conf.type Type of the limits, as in \code{km_ci()}.
#' @param col Fill color of the band.
#'
#' @return \code{NULL}, invisibly.
#'
#' @keywords internal
#' @noRd
km_band_step <- function(km, conf.level, conf.type, col) {
  if (length(km$te) == 0L) return(invisible(NULL))
  ci <- km_ci(km$surv, km$se, conf.level, conf.type)
  xs <- step_extend(km$te, km$tmax)
  lo <- c(1, ci$lower)
  up <- c(1, ci$upper)
  if (length(xs) > length(lo)) {
    lo <- c(lo, lo[length(lo)])
    up <- c(up, up[length(up)])
  }

  # The band runs along the upper step function and back along the lower one
  vl <- step_vertices(xs, lo)
  vu <- step_vertices(xs, up)
  graphics::polygon(c(vu$x, rev(vl$x)), c(vu$y, rev(vl$y)), border = NA,
                    col = col)
  invisible(NULL)
}
