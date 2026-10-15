#' Kaplan-Meier Curves as a Data Frame
#'
#' @description
#' Returns the Kaplan-Meier curves of a \code{kmcurve_fast} object as a data
#' frame, one row per group and distinct observed time, for drawing the curves
#' with other graphics software or for tabulation.
#'
#' @details
#' For each group the rows are the distinct observed times (event and
#' censoring times) in ascending order, preceded by a row at time 0 with
#' survival 1, so that a step function drawn through the rows starts at 1 and
#' ends at the largest observed time. At each time, \code{n.risk} is the
#' number of subjects with an observed time at or after it, \code{n.event} and
#' \code{n.censor} are the numbers of events and censored times at it, and
#' \code{surv} is the Kaplan-Meier estimate just after it. \code{std.err} is
#' the Greenwood standard error on the survival scale, the \code{std.err}
#' reported by \code{summary()} of a \code{survival::survfit()} object, and
#' \code{lower} and \code{upper} are the pointwise confidence limits of the type
#' \code{conf.type}, as drawn by \code{\link{plot.kmcurve_fast}}. At a time at
#' which every subject at risk fails, the estimate is 0 and \code{std.err} and
#' the limits are 0, as in \code{\link{survfit_fast}}, where \code{summary()}
#' reports a \code{NaN} standard error and \code{NA} limits. Apart from the
#' row at time 0, the times, counts, and survival estimates are those reported
#' by \code{summary()} of a \code{survival::survfit()} object with
#' \code{censored = TRUE}.
#'
#' @param x A \code{kmcurve_fast} object from \code{\link{kmcurve_fast}}.
#' @param row.names \code{NULL} or a character vector of row names.
#' @param optional Ignored; the column names are fixed.
#' @param conf.level Confidence level of the pointwise limits. Defaults to
#'   0.95.
#' @param conf.type Type of the pointwise limits: \code{"plain"},
#'   \code{"log"} (default), or \code{"log-log"}, as in
#'   \code{\link{survfit_fast}}.
#' @param ... Further arguments, currently ignored.
#'
#' @return A data frame with the columns \code{group} (the group label, as
#'   printed by \code{\link{print.kmcurve_fast}}), \code{time},
#'   \code{n.risk}, \code{n.event}, \code{n.censor}, \code{surv},
#'   \code{std.err}, \code{lower}, and \code{upper}. In the two-group case the
#'   rows of the control group come first.
#'
#' @seealso \code{\link{kmcurve_fast}}, \code{\link{plot.kmcurve_fast}}
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' t0 <- rexp(n, log(2) / 12)
#' t1 <- rexp(n, log(2) / 18)
#' cens <- runif(2 * n, 0, 30)
#' time  <- pmin(c(t0, t1), cens)
#' event <- as.integer(c(t0, t1) <= cens)
#' group <- rep(c(0, 1), each = n)
#'
#' km <- as.data.frame(kmcurve_fast(time, event, group, control = 0))
#' head(km)
#'
#' @method as.data.frame kmcurve_fast
#' @export
as.data.frame.kmcurve_fast <- function(x, row.names = NULL, optional = FALSE,
                                       conf.level = 0.95, conf.type = "log",
                                       ...) {
  conf.type <- match.arg(conf.type, c("plain", "log", "log-log"))
  if (length(conf.level) != 1L || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be in (0, 1)")
  }

  # Product-limit table of one group at its distinct observed times
  one_group <- function(tt, ee, label) {
    ut       <- sort(unique(tt))
    k        <- length(ut)
    idx      <- match(tt, ut)
    n_event  <- tabulate(idx[ee == 1L], nbins = k)
    n_censor <- tabulate(idx[ee == 0L], nbins = k)
    n_risk   <- rev(cumsum(rev(n_event + n_censor)))
    surv     <- cumprod(1 - n_event / n_risk)
    # Greenwood sum, with no contribution when all subjects at risk fail (as
    # in the curves of kmcurve_fast())
    denom    <- n_risk * (n_risk - n_event)
    gw       <- cumsum(ifelse(denom > 0, n_event / denom, 0))
    std_err  <- surv * sqrt(gw)
    ci       <- km_ci(surv, std_err, conf.level, conf.type)
    data.frame(group    = label,
               time     = c(0, ut),
               n.risk   = c(length(tt), n_risk),
               n.event  = c(0, n_event),
               n.censor = c(0, n_censor),
               surv     = c(1, surv),
               std.err  = c(0, std_err),
               lower    = c(1, ci$lower),
               upper    = c(1, ci$upper),
               stringsAsFactors = FALSE)
  }

  dat <- x$data
  ev  <- as.integer(dat$event)
  out <- if (isTRUE(x$two_group)) {
    is_ctrl <- as.character(dat$group) == as.character(dat$control)
    rbind(one_group(dat$time[is_ctrl], ev[is_ctrl], x$labels[["control"]]),
          one_group(dat$time[!is_ctrl], ev[!is_ctrl], x$labels[["treat"]]))
  } else {
    one_group(dat$time, ev, x$labels[["control"]])
  }
  rownames(out) <- row.names
  out
}
