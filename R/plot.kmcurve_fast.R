#' Plot Method for Kaplan-Meier Curves
#'
#' @description
#' Draws the Kaplan-Meier curves of a \code{kmcurve_fast} object. With
#' \code{rmst = TRUE} the area under each curve up to \code{tau} is shaded to
#' visualize the restricted mean survival time. With \code{hr = TRUE} a second
#' panel shows the smoothed time-varying hazard ratio together with the constant
#' Cox proportional-hazards estimate and, optionally, the average hazard ratio,
#' in the style of stylized nonproportional-hazards figures. For a single-group
#' object only one survival curve is drawn; \code{hr = TRUE} is ignored because a
#' hazard ratio requires two groups.
#'
#' @details
#' The smoothed hazard ratio is obtained by kernel smoothing the Nelson-Aalen
#' increments of each group. For group \code{j} the kernel hazard estimate is
#' the sum over event times of \code{(1 / bw) K((t - t_i) / bw) (d_i / Y_i)},
#' with an Epanechnikov kernel \code{K}, and the hazard ratio is the treatment
#' estimate divided by the control estimate. The bandwidth \code{bw} controls the
#' width of the averaging window and materially affects the result, so it is
#' exposed as an argument. To limit boundary bias the hazard ratio is drawn over
#' an interior time window, and points where the control hazard estimate is
#' effectively zero are dropped.
#'
#' The constant Cox hazard ratio is taken from \code{\link{coxph_fast}}, the
#' average hazard ratio from \code{\link{ahr_fast}}, and the restricted mean
#' survival time difference from \code{\link{rmst_fast}}.
#'
#' @param x A \code{kmcurve_fast} object from \code{\link{kmcurve_fast}}.
#' @param hr Logical. If \code{TRUE}, add a hazard-ratio panel to the right.
#' @param rmst Logical. If \code{TRUE}, shade the area under each curve up to
#'   \code{tau} and annotate the restricted mean survival time difference.
#' @param tau Upper limit used for the restricted mean survival time and the
#'   average hazard ratio. Defaults to the smaller of the two groups' largest
#'   observed times.
#' @param bw Kernel bandwidth for the smoothed hazard ratio. Defaults to a value
#'   derived from the spread of the event times.
#' @param n_grid Number of grid points for the smoothed hazard ratio.
#' @param conf.int Logical. If \code{TRUE}, add pointwise confidence limits to
#'   the survival curves.
#' @param hr_max Upper limit of the hazard-ratio axis. Defaults to a value
#'   derived from the estimated hazard ratio.
#' @param cox_line Logical. If \code{TRUE}, add the constant Cox hazard ratio as
#'   a horizontal line in the hazard-ratio panel.
#' @param ahr_line Logical. If \code{TRUE}, add the average hazard ratio as a
#'   horizontal line in the hazard-ratio panel.
#' @param mark_time Logical. If \code{TRUE}, mark censoring times on the curves.
#' @param conf.level Confidence level for the survival confidence limits.
#' @param conf.type Type of survival confidence limits: one of \code{"plain"},
#'   \code{"log"} (default), or \code{"log-log"}, matching
#'   \code{\link{survfit_fast}}.
#' @param col_control,col_treat,col_hr,col_cox,col_ahr Colors for the control
#'   curve, the treatment curve, the smoothed hazard ratio, the Cox line, and the
#'   average-hazard-ratio line.
#' @param lty_control,lty_treat Line types for the control and treatment curves.
#' @param lty_hr Line type for the smoothed hazard ratio.
#' @param lwd_control,lwd_treat Line widths for the control and treatment curves.
#' @param lwd_hr Line width for the lines in the hazard-ratio panel.
#' @param base_size Relative size of all text.
#' @param legend_pos Keyword giving the legend position in both panels, passed to
#'   \code{\link[graphics]{legend}} (for example \code{"topright"} or
#'   \code{"bottomleft"}).
#' @param xlab,ylab Axis labels for the survival panel.
#' @param main_km,main_hr Panel titles.
#' @param group_labels Optional length-two character vector of legend labels for
#'   the control and treatment groups.
#' @param ... Further arguments, currently ignored.
#'
#' @return The input \code{x}, invisibly.
#'
#' @seealso \code{\link{kmcurve_fast}}, \code{\link{print.kmcurve_fast}}
#'
#' @export
plot.kmcurve_fast <- function(x, hr = FALSE, rmst = FALSE, tau = NULL,
                              bw = NULL, n_grid = 200L, conf.int = TRUE,
                              hr_max = NULL, cox_line = TRUE, ahr_line = FALSE,
                              mark_time = TRUE, conf.level = 0.95,
                              conf.type = "log",
                              col_control = "#1b9e77", col_treat = "#d95f02",
                              col_hr = "#3060c0", col_cox = "black",
                              col_ahr = "grey40",
                              lty_control = 1, lty_treat = 1, lty_hr = 1,
                              lwd_control = 2, lwd_treat = 2, lwd_hr = 2,
                              base_size = 1, legend_pos = "topright",
                              xlab = "Time", ylab = "Survival probability",
                              main_km = "Kaplan-Meier",
                              main_hr = "Hazard ratio",
                              group_labels = NULL, ...) {
  conf.type <- match.arg(conf.type, c("plain", "log", "log-log"))
  two  <- isTRUE(x$two_group)
  km_c <- x$km$control
  km_t <- x$km$treat
  if (!two && isTRUE(hr)) {
    warning("hazard ratio requires two groups; drawing the survival curve only")
    hr <- FALSE
  }
  lab <- if (!is.null(group_labels)) {
    group_labels
  } else if (two) {
    c(paste0("Control (", x$labels[["control"]], ")"),
      paste0("Treatment (", x$labels[["treat"]], ")"))
  } else {
    paste0("Group (", x$labels[["control"]], ")")
  }
  tmax_both <- if (two) min(km_c$tmax, km_t$tmax, na.rm = TRUE) else km_c$tmax
  xmax_km   <- if (two) max(km_c$tmax, km_t$tmax, na.rm = TRUE) else km_c$tmax
  if (is.null(tau)) tau <- tmax_both

  if (isTRUE(hr)) {
    op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 4),
                        cex.lab = base_size, cex.axis = base_size)
  } else {
    op <- graphics::par(mar = c(4, 4, 3, 2),
                        cex.lab = base_size, cex.axis = base_size)
  }
  on.exit(graphics::par(op), add = TRUE)

  ## ---- Kaplan-Meier panel ----
  graphics::plot.default(NA, NA, xlim = c(0, xmax_km), ylim = c(0, 1),
                         xaxs = "i", yaxs = "i",
                         xlab = xlab, ylab = ylab, main = main_km,
                         cex.main = 0.95 * base_size)

  if (isTRUE(rmst)) {
    pc <- step_polygon_to_tau(km_c$te, km_c$surv, tau)
    graphics::polygon(pc$x, pc$y, border = NA,
                      col = grDevices::adjustcolor(col_control, 0.18))
    if (two) {
      pt <- step_polygon_to_tau(km_t$te, km_t$surv, tau)
      graphics::polygon(pt$x, pt$y, border = NA,
                        col = grDevices::adjustcolor(col_treat, 0.18))
    }
    graphics::abline(v = tau, col = "grey60", lty = 3)
  }

  if (isTRUE(conf.int)) {
    km_ci_step(km_c, conf.level, conf.type,
               grDevices::adjustcolor(col_control, 0.5))
    if (two) {
      km_ci_step(km_t, conf.level, conf.type,
                 grDevices::adjustcolor(col_treat, 0.5))
    }
  }

  km_step(km_c, col_control, lty_control, lwd_control)
  if (two) km_step(km_t, col_treat, lty_treat, lwd_treat)

  if (isTRUE(mark_time)) {
    graphics::points(km_c$cens, surv_at(km_c, km_c$cens), pch = 3,
                     col = col_control, cex = 0.6 * base_size)
    if (two) {
      graphics::points(km_t$cens, surv_at(km_t, km_t$cens), pch = 3,
                       col = col_treat, cex = 0.6 * base_size)
    }
  }

  if (isTRUE(rmst)) {
    if (two) {
      rm_res <- tryCatch(
        rmst_fast(x$data$time, x$data$event, group = x$data$group,
                  control = x$data$control, tau = tau),
        error = function(e) NULL)
      if (!is.null(rm_res)) {
        graphics::mtext(sprintf("RMST difference = %.2f at tau = %.1f",
                                unname(rm_res["diff"]), tau),
                        side = 3, line = 0.2, cex = 0.75 * base_size)
      }
    } else {
      rm_res <- tryCatch(
        rmst_fast(x$data$time, x$data$event, tau = tau),
        error = function(e) NULL)
      if (!is.null(rm_res)) {
        graphics::mtext(sprintf("RMST = %.2f at tau = %.1f",
                                unname(rm_res["rmst"]), tau),
                        side = 3, line = 0.2, cex = 0.75 * base_size)
      }
    }
  }

  if (two) {
    legend_fixed(legend_pos, bty = "n", cex = 0.8 * base_size, seg.len = 2.5,
                 legend = lab, col = c(col_control, col_treat),
                 lwd = c(lwd_control, lwd_treat),
                 lty = c(lty_control, lty_treat))
  } else {
    legend_fixed(legend_pos, bty = "n", cex = 0.8 * base_size, seg.len = 2.5,
                 legend = lab, col = col_control,
                 lwd = lwd_control, lty = lty_control)
  }

  ## ---- Hazard-ratio panel ----
  if (isTRUE(hr)) {
    ev_all   <- c(km_c$te, km_t$te)
    hr_upper <- tmax_both
    if (length(ev_all) > 1L) {
      hr_upper <- min(hr_upper,
                      stats::quantile(ev_all, 0.95, names = FALSE,
                                      na.rm = TRUE))
    }
    if (is.null(bw)) {
      bw <- if (length(ev_all) > 1L) 0.15 * diff(range(ev_all)) else hr_upper / 5
      if (!is.finite(bw) || bw <= 0) bw <- hr_upper / 5
    }
    lo <- bw
    hi <- hr_upper - bw
    if (!(hi > lo)) {
      lo <- 0.05 * hr_upper
      hi <- hr_upper
    }
    grid <- seq(lo, hi, length.out = n_grid)
    h_c  <- smooth_hazard(grid, km_c$te, km_c$dLambda, bw)
    h_t  <- smooth_hazard(grid, km_t$te, km_t$dLambda, bw)
    hr_t <- h_t / h_c
    valid <- is.finite(hr_t) & h_c > 1e-8

    if (is.null(hr_max)) {
      hv <- hr_t[valid]
      hr_max <- if (length(hv)) {
        max(1.2, ceiling(stats::quantile(hv, 0.95, names = FALSE) / 0.5) * 0.5)
      } else {
        2
      }
    }
    hr_plot <- hr_t
    hr_plot[!valid | hr_plot > hr_max] <- NA

    graphics::plot.default(grid, hr_plot, type = "l", col = col_hr, lwd = lwd_hr,
                           lty = lty_hr,
                           xlim = c(0, hr_upper), ylim = c(0, hr_max),
                           xaxs = "i", yaxs = "i",
                           xlab = xlab, ylab = "Hazard ratio", main = main_hr,
                           cex.main = 0.95 * base_size)
    graphics::abline(h = 1, col = "grey70", lty = 3)

    leg  <- "Smoothed HR"
    cols <- col_hr
    lwds <- lwd_hr
    ltys <- lty_hr

    if (isTRUE(cox_line)) {
      cx <- tryCatch(
        coxph_fast(x$data$time, x$data$event, x$data$group,
                   control = x$data$control),
        error = function(e) NULL)
      if (!is.null(cx) && is.finite(cx[["exp(coef)"]])) {
        graphics::abline(h = unname(cx[["exp(coef)"]]), col = col_cox,
                         lwd = lwd_hr)
        leg  <- c(leg, "Cox PH")
        cols <- c(cols, col_cox)
        lwds <- c(lwds, lwd_hr)
        ltys <- c(ltys, 1)
      }
    }

    if (isTRUE(ahr_line)) {
      ah <- tryCatch(
        ahr_fast(x$data$time, x$data$event, x$data$group,
                 control = x$data$control, tau = tau),
        error = function(e) NULL)
      if (!is.null(ah) && is.finite(ah$ahr)) {
        graphics::abline(h = ah$ahr, col = col_ahr, lwd = lwd_hr, lty = 2)
        leg  <- c(leg, "Average HR")
        cols <- c(cols, col_ahr)
        lwds <- c(lwds, lwd_hr)
        ltys <- c(ltys, 2)
      }
    }

    t_cross <- hr_cross(grid[valid], hr_t[valid], 1)
    if (!is.null(t_cross)) {
      graphics::segments(t_cross, 0, t_cross, 1, col = col_hr, lty = 3)
      graphics::segments(0, 1, t_cross, 1, col = col_hr, lty = 3)
      graphics::points(t_cross, 1, pch = 19, col = col_hr,
                       cex = 0.9 * base_size)
    }

    legend_fixed(legend_pos, bty = "n", cex = 0.8 * base_size, seg.len = 2.5,
                 legend = leg, col = cols, lwd = lwds, lty = ltys)
  }

  invisible(x)
}

# ------------------------------------------------------------------ #
#  Internal helpers for the plot method
# ------------------------------------------------------------------ #

# Step coordinates for a Kaplan-Meier curve, including the flat segment at 1
# and the flat tail from the last event time to the largest observed time.
km_step <- function(km, col, lty, lwd) {
  xs <- step_extend(km$te, km$tmax)
  ys <- c(1, km$surv)
  if (length(xs) > length(ys)) ys <- c(ys, ys[length(ys)])
  graphics::lines(xs, ys, type = "s", col = col, lty = lty, lwd = lwd)
}

# Pointwise confidence limits drawn as dashed steps, of the requested type.
km_ci_step <- function(km, conf.level, conf.type, col) {
  if (length(km$te) == 0L) return(invisible(NULL))
  ci <- km_ci(km$surv, km$se, conf.level, conf.type)
  xs <- step_extend(km$te, km$tmax)
  lo <- c(1, ci$lower)
  up <- c(1, ci$upper)
  if (length(xs) > length(lo)) {
    lo <- c(lo, lo[length(lo)])
    up <- c(up, up[length(up)])
  }
  graphics::lines(xs, lo, type = "s", col = col, lty = 2)
  graphics::lines(xs, up, type = "s", col = col, lty = 2)
}

# Vectorized Kaplan-Meier confidence limits matching survfit_fast: "plain"
# (survival scale), "log" (default), and "log-log" (complementary log-log).
km_ci <- function(surv, se, conf.level, conf.type) {
  z     <- stats::qnorm(1 - (1 - conf.level) / 2)
  lower <- surv
  upper <- surv
  pos   <- surv > 0 & surv < 1 & is.finite(se)
  if (conf.type == "plain") {
    lower[pos] <- surv[pos] - z * se[pos]
    upper[pos] <- surv[pos] + z * se[pos]
  } else if (conf.type == "log") {
    lower[pos] <- surv[pos] * exp(-z * se[pos] / surv[pos])
    upper[pos] <- surv[pos] * exp( z * se[pos] / surv[pos])
  } else {
    s        <- surv[pos]
    log_s    <- log(s)
    se_ll    <- se[pos] / (s * abs(log_s))
    theta_lo <- log(-log_s) - z * se_ll
    theta_hi <- log(-log_s) + z * se_ll
    lower[pos] <- exp(-exp(theta_hi))
    upper[pos] <- exp(-exp(theta_lo))
  }
  zero <- surv <= 0
  lower[zero] <- 0
  upper[zero] <- 0
  lower <- pmax(0, pmin(1, lower))
  upper <- pmax(0, pmin(1, upper))
  list(lower = lower, upper = upper)
}

# Survival value of a step curve at arbitrary times.
surv_at <- function(km, tt) {
  if (length(tt) == 0L) return(numeric(0))
  if (length(km$te) == 0L) return(rep(1, length(tt)))
  idx <- findInterval(tt, km$te)
  ifelse(idx == 0L, 1, km$surv[pmax(idx, 1L)])
}

# Upper boundary of the area under a step curve up to tau, closed down to zero.
step_polygon_to_tau <- function(te, surv, tau) {
  xs   <- 0
  ys   <- 1
  prev <- 1
  for (i in seq_along(te)) {
    if (te[i] > tau) break
    xs   <- c(xs, te[i], te[i])
    ys   <- c(ys, prev, surv[i])
    prev <- surv[i]
  }
  xs <- c(xs, tau)
  ys <- c(ys, prev)
  list(x = c(xs, tau, 0), y = c(ys, 0, 0))
}

# Epanechnikov kernel.
epan_kernel <- function(u) ifelse(abs(u) <= 1, 0.75 * (1 - u * u), 0)

# Kernel-smoothed hazard from Nelson-Aalen increments.
smooth_hazard <- function(grid, te, dLambda, bw) {
  if (length(te) == 0L) return(rep(NA_real_, length(grid)))
  vapply(grid, function(t0) {
    sum(epan_kernel((t0 - te) / bw) * dLambda) / bw
  }, numeric(1L))
}

# First crossing of a level by a curve, by linear interpolation.
hr_cross <- function(grid, values, level = 1) {
  if (length(grid) < 2L) return(NULL)
  s   <- values - level
  idx <- which(s[-1L] * s[-length(s)] < 0)
  if (length(idx) == 0L) return(NULL)
  i  <- idx[1L]
  t0 <- grid[i]
  t1 <- grid[i + 1L]
  y0 <- s[i]
  y1 <- s[i + 1L]
  t0 - y0 * (t1 - t0) / (y1 - y0)
}
