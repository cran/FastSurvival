#' Plot Method for Scenario Sets
#'
#' @description
#' Draws, for each scenario in a \code{scenario_fast} object, the analytic
#' control and experimental survival curves together with the piecewise-constant
#' hazard ratio on a secondary axis, in the style of the stylized
#' nonproportional-hazards figures used in the design literature. The curves are
#' computed analytically from the piecewise-exponential parameters, so no
#' simulation is involved.
#'
#' @details
#' Each scenario is shown in its own panel arranged on a grid. The hazard ratio
#' is rescaled into the survival plotting region and read on the right-hand axis;
#' a reference line marks a hazard ratio of one. A single-arm scenario (no
#' treatment specification) is drawn as a single survival curve with no hazard
#' ratio.
#'
#' @param x A \code{scenario_fast} object from \code{\link{gen_scenario_fast}}.
#' @param which Optional integer or character vector selecting a subset of
#'   scenarios to draw.
#' @param tmax Upper limit of the time axis. Defaults to a value derived from the
#'   scenarios.
#' @param hr_max Upper limit of the hazard-ratio axis. Defaults to a value
#'   derived from the scenarios.
#' @param n_time Number of points on the time grid.
#' @param mfrow Optional length-two panel layout \code{c(rows, cols)}. Defaults
#'   to a near-square grid.
#' @param col_control,col_treat,col_hr Colors for the control survival curve, the
#'   experimental survival curve, and the hazard ratio.
#' @param lty_control,lty_treat,lty_hr Line types for the control survival curve,
#'   the experimental survival curve, and the hazard ratio.
#' @param lwd_control,lwd_treat,lwd_hr Line widths for the control survival curve,
#'   the experimental survival curve, and the hazard ratio.
#' @param base_size Relative size of all text (titles, axis labels, axis
#'   annotation, and the legend). The default \code{1} reproduces the standard
#'   sizes; values above one enlarge the text and values below one shrink it.
#' @param xlab,ylab Axis labels for the survival panels.
#' @param legend_panel Index (within the drawn subset) of the panel that carries
#'   the legend.
#' @param legend_pos Keyword giving the legend position, passed to
#'   \code{\link[graphics]{legend}} (for example \code{"topright"} or
#'   \code{"bottomleft"}).
#' @param ... Further arguments, currently ignored.
#'
#' @return The input \code{x}, invisibly.
#'
#' @seealso \code{\link{gen_scenario_fast}}, \code{\link{print.scenario_fast}}
#'
#' @export
plot.scenario_fast <- function(x, which = NULL, tmax = NULL, hr_max = NULL,
                               n_time = 481L, mfrow = NULL,
                               col_control = "#1b9e77",
                               col_treat   = "#d95f02",
                               col_hr      = "#3060c0",
                               lty_control = 1,
                               lty_treat   = 1,
                               lty_hr      = 2,
                               lwd_control = 2,
                               lwd_treat   = 2,
                               lwd_hr      = 1.5,
                               base_size   = 1,
                               xlab = "Time",
                               ylab = "Survival probability",
                               legend_panel = 1L,
                               legend_pos   = "topright", ...) {
  scn <- x$scenarios
  if (!is.null(which)) scn <- scn[which]
  n_scn <- length(scn)
  if (n_scn == 0L) stop("no scenarios to plot")

  if (is.null(tmax)) tmax <- default_tmax(scn)
  times <- seq(0, tmax, length.out = n_time)
  evals <- lapply(scn, function(s) scenario_eval(s$args, times))

  if (is.null(hr_max)) {
    hr_all <- unlist(lapply(evals, function(e) e$hr))
    hr_all <- hr_all[is.finite(hr_all)]
    hr_max <- if (length(hr_all)) {
      max(1.2, ceiling(max(hr_all) / 0.2) * 0.2)
    } else {
      1.6
    }
  }

  if (is.null(mfrow)) {
    n_col <- ceiling(sqrt(n_scn))
    n_row <- ceiling(n_scn / n_col)
    mfrow <- c(n_row, n_col)
  }

  op <- graphics::par(mfrow = mfrow, mar = c(4, 4, 3, 4),
                      cex.lab = base_size, cex.axis = base_size)
  on.exit(graphics::par(op), add = TRUE)

  for (k in seq_len(n_scn)) {
    draw_scenario_panel(scn[[k]]$label, evals[[k]], hr_max,
                        col_control, col_treat, col_hr,
                        lty_control, lty_treat, lty_hr,
                        lwd_control, lwd_treat, lwd_hr,
                        base_size, xlab, ylab, legend_pos,
                        add_legend = (k == legend_panel))
  }
  invisible(x)
}

# Draw a single scenario panel.
draw_scenario_panel <- function(label, ev, hr_max,
                                col_control, col_treat, col_hr,
                                lty_control, lty_treat, lty_hr,
                                lwd_control, lwd_treat, lwd_hr,
                                base_size, xlab, ylab, legend_pos,
                                add_legend) {
  graphics::plot.default(ev$times, ev$s_c, type = "l", col = col_control,
                         lwd = lwd_control, lty = lty_control,
                         xaxs = "i", yaxs = "i",
                         ylim = c(0, 1), xlab = xlab,
                         ylab = ylab,
                         main = label, cex.main = 0.95 * base_size)
  if (isTRUE(ev$two)) {
    graphics::lines(ev$times, ev$s_t, col = col_treat,
                    lwd = lwd_treat, lty = lty_treat)
    graphics::par(new = TRUE)
    graphics::plot.default(ev$times, ev$hr / hr_max, type = "l",
                           col = col_hr, lwd = lwd_hr, lty = lty_hr,
                           xaxs = "i", yaxs = "i",
                           ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    graphics::abline(h = 1 / hr_max, col = "grey60", lty = 3)
    hr_ticks <- pretty(c(0, hr_max), n = 4)
    hr_ticks <- hr_ticks[hr_ticks <= hr_max]
    graphics::axis(4, at = hr_ticks / hr_max, labels = hr_ticks)
    graphics::mtext("Hazard ratio", side = 4, line = 2.5, cex = base_size)
  }
  if (add_legend) {
    if (isTRUE(ev$two)) {
      leg  <- c("Control", "Experimental", "Hazard ratio")
      cols <- c(col_control, col_treat, col_hr)
      lwds <- c(lwd_control, lwd_treat, lwd_hr)
      ltys <- c(lty_control, lty_treat, lty_hr)
    } else {
      leg  <- "Survival"
      cols <- col_control
      lwds <- lwd_control
      ltys <- lty_control
    }
    legend_fixed(legend_pos, bty = "n", cex = 0.75 * base_size,
                 seg.len = 2.5, legend = leg, col = cols,
                 lwd = lwds, lty = ltys)
  }
}
