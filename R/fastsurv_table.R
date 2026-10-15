#' Rows of a result of the estimation and testing functions
#'
#' Internal helper for the \code{as.data.frame()}, \code{coef()}, and
#' \code{confint()} methods of the class \code{"fastsurv"}. It returns one row
#' for each estimate and each test that the result reports, in the layout
#' described in \code{\link{as.data.frame.fastsurv}}.
#'
#' @param x A result of class \code{"fastsurv"}.
#'
#' @return A data frame with the columns \code{term}, \code{estimate},
#'   \code{std.error}, \code{conf.low}, \code{conf.high}, \code{statistic},
#'   and \code{p.value}, with the attributes \code{estimates} (a logical
#'   vector marking the rows that are estimates, as opposed to rows of a test
#'   only) and \code{conf.level} (\code{NULL} for a test only).
#'
#' @importFrom stats qnorm pnorm pchisq
#' @keywords internal
#' @noRd
fastsurv_table <- function(x) {
  # Element of a named vector or of a list
  el <- function(nm) unname(unclass(x)[[nm]])

  cl <- if (is.list(x)) x[["conf.level"]] else attr(x, "conf.level")
  q  <- if (is.null(cl)) NA_real_ else qnorm(1 - (1 - cl) / 2)
  side <- if (is.list(x)) x[["side"]] else attr(x, "side")
  if (is.null(side)) side <- 2L

  # Rows with the given columns; a single value is recycled over the rows
  tab <- function(term, estimate = NA, std.error = NA, conf.low = NA,
                  conf.high = NA, statistic = NA, p.value = NA, est = TRUE) {
    k   <- length(term)
    num <- function(v) {
      v <- as.numeric(v)
      if (length(v) == 1L && k > 1L) rep(v, k) else v
    }
    out <- data.frame(term = term, estimate = num(estimate),
                      std.error = num(std.error), conf.low = num(conf.low),
                      conf.high = num(conf.high), statistic = num(statistic),
                      p.value = num(p.value), stringsAsFactors = FALSE)
    attr(out, "estimates")  <- rep(est, k)
    attr(out, "conf.level") <- if (est) cl else NULL
    out
  }

  # Standard error of an estimate whose confidence limits are symmetric on
  # its scale, and the Wald statistic that follows from it
  se_from_limits <- function(lo, hi) (hi - lo) / (2 * q)
  z_from_se <- function(est, se) {
    if (is.finite(est) && is.finite(se) && se > 0) est / se else NA_real_
  }

  if (inherits(x, "survfit_fast")) {
    return(tab("surv", el("surv"), el("std.err"), el("lower"), el("upper")))
  }

  if (inherits(x, "survdiff_fast")) {
    # The value is the signed Z (side = 1) or the chi-square (side = 2); the
    # p-value follows the print method
    stat <- as.numeric(unclass(x))[1L]
    p <- if (is.na(stat)) {
      NA_real_
    } else if (side == 1L) {
      pnorm(stat)
    } else {
      pchisq(stat, df = 1L, lower.tail = FALSE)
    }
    wt <- attr(x, "weight")
    return(tab(if (is.null(wt)) "logrank" else wt, statistic = stat,
               p.value = p, est = FALSE))
  }

  if (inherits(x, "coxph_fast")) {
    # Log hazard ratio, on the scale of its standard error, as coef() and
    # confint() of survival::coxph() report it
    b  <- el("coef")
    se <- el("se(coef)")
    z  <- b / se
    p  <- if (side == 1L) pnorm(z) else 2 * pnorm(-abs(z))
    lo <- log(el(sprintf("lower .%g", cl * 100)))
    hi <- log(el(sprintf("upper .%g", cl * 100)))
    return(tab("coef", b, se, lo, hi, z, p))
  }

  if (inherits(x, "rmst_fast")) {
    if (!("diff" %in% names(x))) {
      return(tab("rmst", el("rmst"), el("std.err"), el("lower"), el("upper")))
    }
    # The ratio is reported on the log scale, on which its interval is
    # symmetric and its test is computed
    lo_r <- log(el("ratio.lower"))
    hi_r <- log(el("ratio.upper"))
    return(tab(c("rmst.ctrl", "rmst.trt", "diff", "log.ratio"),
               c(el("rmst.ctrl"), el("rmst.trt"), el("diff"),
                 log(el("ratio"))),
               c(NA, NA, el("se.diff"), se_from_limits(lo_r, hi_r)),
               c(NA, NA, el("diff.lower"), lo_r),
               c(NA, NA, el("diff.upper"), hi_r),
               c(NA, NA, el("z.diff"), el("z.ratio")),
               c(NA, NA, el("p.diff"), el("p.ratio"))))
  }

  if (inherits(x, "wmst_fast") || inherits(x, "medsurv_fast")) {
    stem <- if (inherits(x, "wmst_fast")) "wmst" else "median"
    if (!isTRUE(attr(x, "two.group"))) {
      return(tab(stem, el(stem), el("se"), el("lower"), el("upper")))
    }
    return(tab(c(paste0(stem, c(".control", ".treatment")), "diff"),
               c(el(paste0(stem, ".control")), el(paste0(stem, ".treatment")),
                 el("diff")),
               c(el("se.control"), el("se.treatment"), el("se.diff")),
               c(el("lower.control"), el("lower.treatment"),
                 el("lower.diff")),
               c(el("upper.control"), el("upper.treatment"),
                 el("upper.diff")),
               c(NA, NA, el("z")),
               c(NA, NA, el("p"))))
  }

  if (inherits(x, "milestone_fast")) {
    s  <- x[["surv"]]
    se <- x[["std.err"]]
    # The Wald interval of the difference uses the unpooled Greenwood standard
    # error; the "loglog" and "mover" intervals are not of the form
    # estimate +/- q * se, so no standard error is given for them
    se_d <- if (identical(x[["method"]], "wald")) {
      sqrt(se[["control"]]^2 + se[["treatment"]]^2)
    } else {
      NA_real_
    }
    return(tab(c("surv.control", "surv.treatment", "diff"),
               c(s[["control"]], s[["treatment"]], x[["diff"]]),
               c(se[["control"]], se[["treatment"]], se_d),
               c(x[["surv.lower"]][["control"]],
                 x[["surv.lower"]][["treatment"]], x[["diff.lower"]]),
               c(x[["surv.upper"]][["control"]],
                 x[["surv.upper"]][["treatment"]], x[["diff.upper"]]),
               c(NA, NA, x[["statistic"]]),
               c(NA, NA, x[["p.value"]])))
  }

  if (inherits(x, "maxcombo_fast")) {
    z <- attr(x, "z")
    return(tab(c("maxcombo", names(z)),
               statistic = c(el("statistic"), unname(z)),
               p.value = c(el("p.value"), rep(NA_real_, length(z))),
               est = FALSE))
  }

  if (inherits(x, "rmw_fast")) {
    z <- attr(x, "z")
    return(tab(c("rmw", "logrank", "mwlrt"),
               statistic = c(el("statistic"), z[["logrank"]], z[["mwlrt"]]),
               p.value = c(el("p.value"), NA, NA), est = FALSE))
  }

  if (inherits(x, "wkm_fast")) {
    return(tab("wdiff", el("wdiff"), el("se"), el("lower"), el("upper"),
               el("z"), el("p")))
  }

  if (inherits(x, "ahsw_fast")) {
    # The ratio is reported on the log scale; the standard errors of the log
    # ratio and of the difference, which ahsw_fast() does not store, follow
    # from the confidence limits, which are symmetric on these scales
    lo_r <- log(el("rah.lower"))
    hi_r <- log(el("rah.upper"))
    se_r <- se_from_limits(lo_r, hi_r)
    se_d <- se_from_limits(el("dah.lower"), el("dah.upper"))
    return(tab(c("ah.ctrl", "ah.trt", "log.rah", "dah"),
               c(el("ah.ctrl"), el("ah.trt"), log(el("rah")), el("dah")),
               c(NA, NA, se_r, se_d),
               c(NA, NA, lo_r, el("dah.lower")),
               c(NA, NA, hi_r, el("dah.upper")),
               c(NA, NA, z_from_se(log(el("rah")), se_r),
                 z_from_se(el("dah"), se_d)),
               c(NA, NA, el("p.rah"), el("p.dah"))))
  }

  if (inherits(x, "ahr_fast")) {
    # The log average hazard ratio with its test on the log scale, and the
    # treatment group's share of the hazard with the primary test
    return(tab(c("log.ahr", "theta"),
               c(x[["log.ahr"]], unname(x[["theta"]][2L])),
               c(x[["se.loghr"]], x[["se.theta"]]),
               c(log(x[["lower"]]), NA),
               c(log(x[["upper"]]), NA),
               c(x[["z.loghr"]], x[["z"]]),
               c(x[["p.value.loghr"]], x[["p.value"]])))
  }

  stop("no table is defined for an object of class '", class(x)[1L], "'")
}
