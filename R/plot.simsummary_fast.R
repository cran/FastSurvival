#' Plot Method for Sequential Analysis Summaries
#'
#' @description
#' Draws the stopping probabilities of a \code{\link{simsummary_fast}} result
#' by look. With \code{what = "cumulative"} (default) one panel shows, for
#' every block (arm or population), the cumulative probability of stopping
#' for efficacy, which at the final look is the rejection rate, and, when
#' futility boundaries were applied, the cumulative probability of stopping
#' for futility. With \code{what = "stage"} each block has a panel of stacked
#' bars of the probabilities of stopping for efficacy and for futility at each
#' look.
#'
#' @details
#' The plot uses the look rows of each block, which must keep all the looks of
#' the design in their original order and equal the rows of the stored result,
#' as in the print method (see \code{\link{print.simsummary_fast}}): a
#' selection of columns or of some of the looks, rows in a different order, or
#' rows of another result combined by \code{rbind()} give an error. The overall
#' rows are not used, so they are not compared. A result created by
#' FastSurvival 1.1.0, which does not store the looks and rows, gives an error
#' that asks for \code{simsummary_fast()} to be run again. The cumulative
#' probability of stopping for futility is the cumulative sum of
#' \code{prob.stop.futility} over the looks.
#'
#' The looks are placed on the horizontal axis by their index
#' (\code{xaxis = "look"}), by the expected number of events at the look
#' (\code{"events"}, the column \code{n.event.mean}), or by the expected
#' calendar time of the look (\code{"time"}, the column \code{cutoff.mean}).
#'
#' @param x An object of class \code{"simsummary_fast"} from
#'   \code{\link{simsummary_fast}}.
#' @param what \code{"cumulative"} (default) for the cumulative stopping
#'   probabilities of all blocks in one panel, or \code{"stage"} for the
#'   stopping probabilities at each look, one panel per block.
#' @param xaxis The horizontal axis: \code{"look"} (default), \code{"events"},
#'   or \code{"time"}. See Details.
#' @param col Colors of the blocks, recycled. Defaults to a qualitative
#'   palette.
#' @param ylim Limits of the vertical axis. Defaults to the range from 0 to
#'   slightly above the largest probability drawn.
#' @param legend_pos Keyword giving the legend position, passed to
#'   \code{\link[graphics]{legend}}.
#' @param xlab,ylab,main Axis labels and title. Defaults depend on \code{what}
#'   and \code{xaxis}; with \code{what = "stage"} and several blocks, the
#'   panels are titled by their blocks.
#' @param ... Further arguments, currently ignored.
#'
#' @return The input \code{x}, invisibly.
#'
#' @seealso \code{\link{simsummary_fast}}, \code{\link{print.simsummary_fast}}.
#'
#' @examples
#' df <- simdata_fast(
#'   nsim     = 200,
#'   n        = c(150, 150),
#'   a.time   = c(0, 12),
#'   a.rate   = 300 / 12,
#'   e.hazard = list(0.05, 0.035),
#'   seed     = 1
#' )
#' res <- analysis_fast(df, control = 1, event.looks = c(60, 105, 150),
#'                      stat = c("logrank", "coxph"), side = 1)
#' oc <- simsummary_fast(res,
#'                       eff.col = "logrank.z",
#'                       efficacy = c(NA, -2.96, -1.97),
#'                       fut.col = "cox.coef",
#'                       futility = c(log(1.2), NA, NA),
#'                       direction = "lower")
#' plot(oc)
#' plot(oc, what = "stage", xaxis = "events")
#'
#' @method plot simsummary_fast
#' @export
plot.simsummary_fast <- function(x, what = c("cumulative", "stage"),
                                 xaxis = c("look", "events", "time"),
                                 col = NULL, ylim = NULL,
                                 legend_pos = "topleft",
                                 xlab = NULL, ylab = NULL, main = NULL, ...) {
  what  <- match.arg(what)
  xaxis <- match.arg(xaxis)

  # A result of FastSurvival 1.1.0 stores the boundary settings without the
  # looks and rows that the blocks are checked against
  bd <- attr(x, "boundary")
  if (!is.null(bd) && is.null(bd$looks)) {
    stop("this result was created by an earlier version of FastSurvival; ",
         "call simsummary_fast() again to plot it")
  }
  blocks <- simsummary_blocks(x)
  if (is.null(blocks)) {
    stop("plot() needs a simsummary_fast() result whose blocks keep all ",
         "their looks in the original order")
  }
  has_fut <- identical(bd$mode, "z") && !is.null(bd$futility)

  xcol <- switch(xaxis, look = NULL, events = "n.event.mean",
                 time = "cutoff.mean")
  if (!is.null(xcol) && !(xcol %in% names(blocks[[1L]]))) {
    stop("xaxis = \"", xaxis, "\" needs the column '", xcol, "', which the ",
         "result does not have")
  }
  if (!is.null(xcol) &&
      !any(is.finite(unlist(lapply(blocks, function(b) b[[xcol]]))))) {
    stop("xaxis = \"", xaxis, "\" needs finite values in the column '", xcol,
         "'")
  }
  if (is.null(xlab)) {
    xlab <- switch(xaxis, look = "Look", events = "Expected number of events",
                   time = "Expected calendar time")
  }

  nb   <- length(blocks)
  labs <- vapply(blocks, function(b) attr(b, "label"), character(1L))
  if (is.null(col)) {
    col <- if (nb == 1L) "#3060c0" else grDevices::hcl.colors(nb, "Dark 3")
  }
  col <- rep_len(col, nb)

  # Abscissae of the looks of a block, and the cumulative probabilities
  x_of <- function(b) {
    if (is.null(xcol)) seq_len(nrow(b)) else b[[xcol]]
  }
  cum_eff <- lapply(blocks, function(b) b$cum.reject)
  cum_fut <- lapply(blocks, function(b) cumsum(b$prob.stop.futility))

  top_of <- function(v) {
    m <- suppressWarnings(max(unlist(v), na.rm = TRUE))
    if (!is.finite(m) || m <= 0) 1 else min(1, 1.08 * m)
  }

  if (what == "cumulative") {
    if (is.null(ylim)) {
      ylim <- c(0, top_of(c(cum_eff, if (has_fut) cum_fut)))
    }
    if (is.null(ylab)) ylab <- "Cumulative probability of stopping"
    if (is.null(main)) main <- "Cumulative stopping probability"
    xs   <- unlist(lapply(blocks, x_of))
    xlim <- range(xs[is.finite(xs)])
    if (!(diff(xlim) > 0)) xlim <- xlim + c(-0.5, 0.5)

    graphics::plot.default(NA, NA, xlim = xlim, ylim = ylim, xaxt = "n",
                           xlab = xlab, ylab = ylab, main = main)
    if (is.null(xcol)) {
      graphics::axis(1L, at = seq_along(bd$looks), labels = bd$looks)
    } else {
      graphics::axis(1L)
    }
    for (i in seq_len(nb)) {
      xi <- x_of(blocks[[i]])
      graphics::lines(xi, cum_eff[[i]], type = "b", pch = 19, lty = 1,
                      col = col[i], lwd = 2)
      if (has_fut) {
        graphics::lines(xi, cum_fut[[i]], type = "b", pch = 1, lty = 2,
                        col = col[i], lwd = 2)
      }
    }

    # Legend: blocks by color, and efficacy and futility by line type
    leg <- character(0)
    lc  <- character(0)
    lt  <- numeric(0)
    lp  <- numeric(0)
    if (nb > 1L) {
      leg <- labs
      lc  <- col
      lt  <- rep(1, nb)
      lp  <- rep(NA, nb)
    }
    blk_col <- if (nb > 1L) "black" else col[1L]
    leg <- c(leg, "Efficacy")
    lc  <- c(lc, blk_col)
    lt  <- c(lt, 1)
    lp  <- c(lp, 19)
    if (has_fut) {
      leg <- c(leg, "Futility")
      lc  <- c(lc, blk_col)
      lt  <- c(lt, 2)
      lp  <- c(lp, 1)
    }
    graphics::legend(legend_pos, legend = leg, col = lc, lty = lt, pch = lp,
                     lwd = 2, bty = "n", cex = 0.85)
    return(invisible(x))
  }

  # what == "stage": stacked bars of the stopping probabilities at each look
  if (is.null(ylab)) ylab <- "Probability of stopping"
  if (is.null(ylim)) {
    ylim <- c(0, top_of(lapply(blocks, function(b) {
      b$prob.stop.efficacy + b$prob.stop.futility
    })))
  }
  if (nb > 1L) {
    op <- graphics::par(mfrow = grDevices::n2mfrow(nb))
    on.exit(graphics::par(op), add = TRUE)
  }
  for (i in seq_len(nb)) {
    b <- blocks[[i]]
    h <- if (has_fut) {
      rbind(b$prob.stop.efficacy, b$prob.stop.futility)
    } else {
      matrix(b$prob.stop.efficacy, nrow = 1L)
    }
    nms <- if (is.null(xcol)) {
      bd$looks
    } else {
      format(round(b[[xcol]], 1), trim = TRUE)
    }
    ttl <- if (!is.null(main)) {
      main
    } else if (nb > 1L) {
      labs[i]
    } else {
      "Stopping probability by look"
    }
    fills <- c(col[i], grDevices::adjustcolor(col[i], 0.35))[seq_len(nrow(h))]
    graphics::barplot(h, names.arg = nms, col = fills, border = NA,
                      ylim = ylim, xlab = xlab, ylab = ylab, main = ttl)
    if (i == 1L) {
      graphics::legend(legend_pos, legend = c("Efficacy", "Futility")[
                         seq_len(nrow(h))],
                       fill = fills, border = NA, bty = "n", cex = 0.85)
    }
  }
  invisible(x)
}
