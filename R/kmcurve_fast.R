#' Build Kaplan-Meier Curves for One Trial Realization
#'
#' @description
#' Builds the Kaplan-Meier step functions for a single realized data set (for
#' example one replicate from \code{\link{simdata_fast}}) and returns a
#' \code{kmcurve_fast} object that can be drawn with
#' \code{\link{plot.kmcurve_fast}} and summarized with
#' \code{\link{print.kmcurve_fast}}. Both the two-group case and the single-group
#' case are supported. This is an analysis-stage and plotting helper, written in
#' plain R for use on a single data set, and is not intended to be called inside
#' a simulation loop. The fast point estimator \code{\link{survfit_fast}} remains
#' the tool for repeated evaluation inside loops.
#'
#' @details
#' For each group the object stores the distinct event times, the number of
#' events and the number at risk at those times, the Kaplan-Meier survival
#' estimate, the Greenwood standard error, the Nelson-Aalen increments
#' \code{d_i / Y_i} (used by the plot method to estimate a smoothed time-varying
#' hazard ratio), and the censoring times. The pooled data are also retained so
#' that the plot and print methods can reuse \code{\link{coxph_fast}},
#' \code{\link{rmst_fast}}, and \code{\link{ahr_fast}} for the constant Cox
#' hazard ratio, the restricted mean survival time, and the average hazard
#' ratio. When \code{group} is \code{NULL} or has a single distinct value, a
#' single-group object is built and the plot method draws one survival curve
#' with no hazard ratio.
#'
#' @param time A numeric vector of follow-up times for all subjects.
#' @param event An integer or numeric vector of event indicators (1 = event,
#'   0 = censored), aligned with \code{time}.
#' @param group A vector identifying the groups, aligned with \code{time}, with
#'   one or two distinct values. If \code{NULL} (default), all subjects are
#'   treated as a single group.
#' @param control The value of \code{group} that denotes the reference (control)
#'   group; the other value is treated as the treatment group. Required in the
#'   two-group case and ignored in the single-group case.
#'
#' @return An object of class \code{kmcurve_fast}: a list with a logical
#'   \code{two_group} flag, the group values and their labels, the per-group
#'   Kaplan-Meier summaries in element \code{km}, the pooled \code{data}, and the
#'   per-group sample sizes \code{n} and event counts \code{events}.
#'
#' @examples
#' set.seed(1)
#' n <- 150
#' t0 <- rexp(n, log(2) / 12)
#' t1 <- rexp(n, log(2) / 18)
#' cens <- runif(2 * n, 0, 30)
#' time  <- pmin(c(t0, t1), cens)
#' event <- as.integer(c(t0, t1) <= cens)
#' group <- rep(c(0, 1), each = n)
#'
#' fit <- kmcurve_fast(time, event, group, control = 0)
#' print(fit)
#' plot(fit)
#' plot(fit, hr = TRUE, rmst = TRUE, tau = 18)
#'
#' # Single-group object: pass group = NULL (or a single-level group)
#' fit1 <- kmcurve_fast(time[group == 0], event[group == 0])
#' plot(fit1, rmst = TRUE, tau = 18)
#'
#' @seealso \code{\link{plot.kmcurve_fast}}, \code{\link{print.kmcurve_fast}},
#'   \code{\link{survfit_fast}}, \code{\link{coxph_fast}},
#'   \code{\link{rmst_fast}}, \code{\link{ahr_fast}}
#'
#' @export
kmcurve_fast <- function(time, event, group = NULL, control = NULL) {
  n <- length(time)
  event <- as.integer(event)
  if (length(event) != n) {
    stop("'time' and 'event' must have the same length")
  }

  # Determine the number of groups
  if (is.null(group)) {
    ulev  <- NA
    n_grp <- 1L
  } else {
    if (length(group) != n) {
      stop("'group' must have the same length as 'time'")
    }
    if (is.factor(group)) group <- as.character(group)
    ulev  <- sort(unique(group))
    n_grp <- length(ulev)
  }
  if (n_grp > 2L) {
    stop("'group' must have at most two distinct values")
  }

  # ---- Single-group mode ----
  if (n_grp <= 1L) {
    km_s      <- km_one_group(time, event)
    grp_label <- if (is.null(group)) "All" else as.character(ulev)
    return(structure(
      list(
        two_group = FALSE,
        control   = if (is.null(group)) NA else ulev,
        treat     = NA,
        labels    = c(control = grp_label, treat = NA_character_),
        km        = list(control = km_s, treat = NULL),
        data      = list(time = time, event = event, group = group,
                         control = if (is.null(group)) NULL else ulev),
        n         = c(control = km_s$n),
        events    = c(control = km_s$n_event)
      ),
      class = "kmcurve_fast"
    ))
  }

  # ---- Two-group mode ----
  if (is.null(control) || length(control) != 1L ||
      !(as.character(control) %in% as.character(ulev))) {
    stop("'control' must be one of the two values in 'group'")
  }
  ctrl_val <- ulev[as.character(ulev) == as.character(control)]
  trt_val  <- ulev[as.character(ulev) != as.character(control)]

  is_ctrl <- as.character(group) == as.character(ctrl_val)
  km_c <- km_one_group(time[is_ctrl],  event[is_ctrl])
  km_t <- km_one_group(time[!is_ctrl], event[!is_ctrl])

  structure(
    list(
      two_group = TRUE,
      control = ctrl_val,
      treat   = trt_val,
      labels  = c(control = as.character(ctrl_val),
                  treat   = as.character(trt_val)),
      km      = list(control = km_c, treat = km_t),
      data    = list(time = time, event = event, group = group,
                     control = ctrl_val),
      n       = c(control = km_c$n, treat = km_t$n),
      events  = c(control = km_c$n_event, treat = km_t$n_event)
    ),
    class = "kmcurve_fast"
  )
}

# ------------------------------------------------------------------ #
#  Internal helper: single-group Kaplan-Meier and Nelson-Aalen pieces
# ------------------------------------------------------------------ #

# Product-limit estimate with Greenwood standard error and Nelson-Aalen
# increments for one group, computed in plain R.
km_one_group <- function(t_in, e_in) {
  n   <- length(t_in)
  ord <- order(t_in)
  tt  <- t_in[ord]
  ee  <- as.integer(e_in[ord])
  ev_t <- tt[ee == 1L]
  cens <- tt[ee == 0L]

  if (length(ev_t) == 0L) {
    return(list(te = numeric(0), d = numeric(0), Y = numeric(0),
                dLambda = numeric(0), surv = numeric(0), se = numeric(0),
                cens = cens, n = n, n_event = 0L,
                tmax = if (n > 0L) max(tt) else NA_real_))
  }

  te <- sort(unique(ev_t))
  k  <- length(te)
  d  <- as.numeric(tabulate(match(ev_t, te), nbins = k))
  Y  <- vapply(te, function(u) sum(tt >= u), numeric(1L))

  dLambda <- d / Y
  surv    <- cumprod(1 - dLambda)
  denom   <- Y * (Y - d)
  gw_term <- ifelse(denom > 0, d / denom, 0)
  se      <- surv * sqrt(cumsum(gw_term))

  list(te = te, d = d, Y = Y, dLambda = dLambda, surv = surv, se = se,
       cens = cens, n = n, n_event = length(ev_t), tmax = max(tt))
}

# Internal helper: place a legend at a keyword position computed explicitly from
# par("usr"). Keyword positions passed straight to legend() can misbehave under
# some multi-panel (mfrow) device states, so the corner coordinates and
# justification are derived directly and passed as x, y, xjust, yjust.
legend_fixed <- function(pos, inset = 0.04, ...) {
  usr <- graphics::par("usr")
  dx  <- inset * (usr[2L] - usr[1L])
  dy  <- inset * (usr[4L] - usr[3L])
  xl  <- usr[1L] + dx
  xr  <- usr[2L] - dx
  xc  <- (usr[1L] + usr[2L]) / 2
  yb  <- usr[3L] + dy
  yt  <- usr[4L] - dy
  yc  <- (usr[3L] + usr[4L]) / 2
  spec <- switch(
    pos,
    "topright"    = list(x = xr, y = yt, xj = 1,   yj = 1),
    "topleft"     = list(x = xl, y = yt, xj = 0,   yj = 1),
    "bottomright" = list(x = xr, y = yb, xj = 1,   yj = 0),
    "bottomleft"  = list(x = xl, y = yb, xj = 0,   yj = 0),
    "top"         = list(x = xc, y = yt, xj = 0.5, yj = 1),
    "bottom"      = list(x = xc, y = yb, xj = 0.5, yj = 0),
    "left"        = list(x = xl, y = yc, xj = 0,   yj = 0.5),
    "right"       = list(x = xr, y = yc, xj = 1,   yj = 0.5),
    "center"      = list(x = xc, y = yc, xj = 0.5, yj = 0.5),
    NULL)
  if (is.null(spec)) {
    graphics::legend(pos, ...)
  } else {
    graphics::legend(x = spec$x, y = spec$y, xjust = spec$xj, yjust = spec$yj,
                     ...)
  }
}
