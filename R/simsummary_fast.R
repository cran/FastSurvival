#' Summarize Operating Characteristics from Sequential Analysis Output
#'
#' @description
#' Aggregates the per-simulation, per-look output of \code{\link{analysis_fast}}
#' into operating characteristics: the rejection rate, the futility-stopping
#' rate, the distribution of the stopping look, and the expected analysis timing
#' (events and calendar time at stopping). A group-sequential design is summarized
#' by applying the supplied per-look boundaries in sequence, and a fixed design is
#' the single-look (\code{K = 1}) degenerate case of the same logic. This function
#' consumes boundaries computed elsewhere (for example by gsDesign or rpact) and
#' does not compute or spend alpha itself; it estimates the stopping probabilities
#' by Monte Carlo over the simulated trials rather than by the analytic numerical
#' integration used by those packages, so the two agree only up to Monte Carlo
#' error and converge as the number of simulations grows.
#'
#' @details
#' The boundaries are supplied per look, one value for each distinct value of the
#' \code{look} column in \code{data}. Two boundary modes are available and exactly
#' one must be used.
#'
#' In the Z mode the efficacy boundary is compared with the column named by
#' \code{eff.col} and the futility boundary with the column named by
#' \code{fut.col}. The two columns may differ, which lets the efficacy and
#' futility rules live on different scales. For example a beta-spending efficacy
#' boundary can be applied to the standardized statistic \code{logrank.z} while a
#' futility boundary expressed as a log hazard ratio is applied to
#' \code{cox.coef}. The crossing direction is set by \code{direction}: with
#' \code{"lower"} (the natural sign of \code{logrank.z}, \code{cox.z}, and
#' \code{cox.coef}, where treatment benefit is negative) the efficacy boundary is
#' crossed when the efficacy statistic is at or below \code{efficacy} and the
#' futility boundary is crossed when the futility statistic is at or above
#' \code{futility}; with \code{"upper"} the inequalities are reversed. Either
#' boundary vector may contain \code{NA} at some looks to omit that rule there, so
#' efficacy-only and futility-only looks are expressed by placing \code{NA} in the
#' other vector. A look with \code{NA} on both rules can never stop the trial.
#'
#' In the p mode the p-value in the column named by \code{p.col} is compared with
#' the per-look nominal level \code{alpha}, rejecting when the p-value is at or
#' below the level. Futility is not used in the p mode, and \code{alpha} may
#' contain \code{NA} to omit the efficacy test at a look.
#'
#' For each simulated trial the looks are examined in order. The trial stops at
#' the first look whose statistic crosses a boundary. Crossing the efficacy
#' boundary is a rejection of the null hypothesis; crossing the futility boundary
#' (Z mode only) is a stop without rejection. When both are crossed at the same
#' look the efficacy stop takes precedence. A look whose relevant statistic is
#' \code{NA}, or whose boundary is \code{NA}, triggers no crossing of that rule
#' and the trial continues. A trial that reaches the final look without crossing
#' the efficacy boundary does not reject.
#'
#' Each look's \code{prob.stop.efficacy} is the marginal probability of stopping
#' for efficacy for the first time at that look, that is the probability of not
#' stopping at any earlier look and crossing the efficacy boundary at this one.
#' This is the stage-wise rejection contribution of a group-sequential design; the
#' cumulative sum \code{cum.reject} is the cumulative power up to and including
#' that look, and the total over all looks is \code{rejection.rate}. These are the
#' same quantities that gsDesign and rpact report, estimated here by simulation.
#'
#' The rejection rate is the type I error under a null data-generating truth and
#' the power under an alternative truth, but because the function does not know the
#' truth used to generate \code{data} it is reported neutrally as the rejection
#' rate and its interpretation is left to the user.
#'
#' When \code{data} carries a \code{population} column (the long form produced by
#' \code{\link{analysis_fast}} with \code{by.subgroup = TRUE}), the same
#' boundaries are applied within each population and the output has one block of
#' rows per population.
#'
#' @param data A data frame from \code{\link{analysis_fast}}, containing at least
#'   \code{sim} and \code{look}, the statistic columns named by the boundary
#'   arguments, and, when available, \code{n.event} and \code{cutoff} for the
#'   analysis-timing summaries.
#' @param eff.col A single character naming the statistic column compared with
#'   \code{efficacy}. Required with the Z mode.
#' @param efficacy A numeric vector of efficacy boundaries, one per look, on the
#'   scale of \code{eff.col}. Entries may be \code{NA} to omit the efficacy test
#'   at that look. Required with the Z mode.
#' @param fut.col A single character naming the statistic column compared with
#'   \code{futility}. Defaults to \code{eff.col}. Used only in the Z mode and only
#'   when \code{futility} is supplied; it may differ from \code{eff.col} so that
#'   futility is judged on a different scale (for example a log hazard ratio).
#' @param futility A numeric vector of futility boundaries, one per look, on the
#'   scale of \code{fut.col}, or \code{NULL} for no futility stopping. Entries may
#'   be \code{NA} to omit a futility boundary at that look. Used only in the Z
#'   mode.
#' @param direction A single string, either \code{"lower"} (default) or
#'   \code{"upper"}, giving the direction in which the boundaries are crossed. See
#'   Details.
#' @param p.col A single character naming the p-value column compared with
#'   \code{alpha}. Selects the p mode; use either the Z mode (\code{eff.col} with
#'   \code{efficacy}) or the p mode (\code{p.col} with \code{alpha}), not both.
#' @param alpha A numeric vector of per-look nominal significance levels, one per
#'   look. Entries may be \code{NA} to omit the test at a look. Required with
#'   \code{p.col}.
#'
#' @return An object of class \code{"simsummary_fast"}: a data frame with one row
#'   per population and look plus an \code{overall} summary row appended after each
#'   population's looks. The columns are \code{population}, \code{look} (the look
#'   index, or \code{"overall"} on the summary row), optionally \code{look.value},
#'   \code{n.enrolled.mean} and \code{n.event.mean} (the mean enrolled and event
#'   counts at that look, or at the stopping look on the summary row),
#'   \code{n.dropout.mean} and \code{n.pipeline.mean} (the mean dropout count and
#'   pipeline count \code{n.enrolled - n.event - n.dropout}, when those columns
#'   are present in \code{data}), \code{cutoff.mean} (the mean calendar time,
#'   likewise), \code{prob.stop.efficacy},
#'   \code{prob.stop.futility}, \code{prob.stop.any}, and \code{cum.reject}. On the
#'   \code{overall} row \code{prob.stop.efficacy} is the total rejection rate,
#'   \code{prob.stop.futility} the total futility rate, \code{prob.stop.any} their
#'   sum, and \code{cum.reject} again the total rejection rate; its timing columns
#'   are the expected counts and calendar time at the stopping look. The number of
#'   simulations is stored in the attribute \code{nsim} and the boundary settings
#'   in the attribute \code{boundary}.
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
#'
#' res <- analysis_fast(df, control = 1, event.looks = c(60, 105, 150),
#'                      stat = c("logrank", "coxph"), side = 1)
#'
#' # Efficacy on the standardized log-rank Z, futility on the log hazard ratio,
#' # with a futility-only first look and efficacy-only later looks
#' simsummary_fast(res,
#'                 eff.col = "logrank.z",
#'                 efficacy = c(NA, -2.96, -1.97),
#'                 fut.col = "cox.coef",
#'                 futility = c(log(1.2), NA, NA),
#'                 direction = "lower")
#'
#' # p-value boundaries instead
#' simsummary_fast(res, p.col = "logrank.p",
#'                 alpha = c(0.0006, 0.0151, 0.0245))
#'
#' @seealso \code{\link{analysis_fast}}, \code{\link{simdata_fast}},
#'   \code{\link{print.simsummary_fast}}.
#'
#' @export
simsummary_fast <- function(data,
                            eff.col = NULL, efficacy = NULL,
                            fut.col = eff.col, futility = NULL,
                            direction = c("lower", "upper"),
                            p.col = NULL, alpha = NULL) {

  direction <- match.arg(direction)

  if (!is.data.frame(data) || !all(c("sim", "look") %in% names(data))) {
    stop("'data' must be a data frame from analysis_fast with 'sim' and 'look' ",
         "columns")
  }

  use_z <- !is.null(eff.col)
  use_p <- !is.null(p.col)
  if (use_z == use_p) {
    stop("Supply exactly one of 'eff.col' (with 'efficacy') and 'p.col' ",
         "(with 'alpha')")
  }

  looks  <- sort(unique(data$look))
  n_look <- length(looks)

  if (use_z) {
    if (length(eff.col) != 1L || !eff.col %in% names(data)) {
      stop("'eff.col' must name a single column of 'data'")
    }
    if (is.null(efficacy) || length(efficacy) != n_look) {
      stop("'efficacy' must have one boundary per look (", n_look, ")")
    }
    use_fut <- !is.null(futility)
    if (use_fut) {
      if (length(fut.col) != 1L || !fut.col %in% names(data)) {
        stop("'fut.col' must name a single column of 'data'")
      }
      if (length(futility) != n_look) {
        stop("'futility' must have one boundary per look (", n_look,
             ") or be NULL")
      }
    }
  } else {
    if (length(p.col) != 1L || !p.col %in% names(data)) {
      stop("'p.col' must name a single column of 'data'")
    }
    if (is.null(alpha) || length(alpha) != n_look) {
      stop("'alpha' must have one nominal level per look (", n_look, ")")
    }
    if (!is.null(futility)) {
      stop("'futility' is supported only with the Z mode ('eff.col')")
    }
    use_fut <- FALSE
  }

  has_enrol   <- "n.enrolled" %in% names(data)
  has_event   <- "n.event" %in% names(data)
  has_dropout <- "n.dropout" %in% names(data)
  has_pipeline <- "n.pipeline" %in% names(data)
  has_cutoff  <- "cutoff" %in% names(data)
  has_lookval <- "look.value" %in% names(data)
  has_pop     <- "population" %in% names(data)
  pops        <- if (has_pop) unique(data$population) else "overall"

  # Boundary description stored as an attribute for the print method
  bdesc <- if (use_z) {
    list(mode = "z", eff.col = eff.col, efficacy = efficacy,
         fut.col = if (use_fut) fut.col else NA_character_,
         futility = futility, direction = direction)
  } else {
    list(mode = "p", p.col = p.col, alpha = alpha)
  }

  # ------------------------------------------------------------------ #
  #  Summarize one population's rows (all sims, all looks)
  # ------------------------------------------------------------------ #
  summarize_one <- function(df, pop_label) {
    sims <- sort(unique(df$sim))
    nsim <- length(sims)
    si   <- match(df$sim, sims)
    li   <- match(df$look, looks)

    fill_mat <- function(col) {
      m <- matrix(NA_real_, nsim, n_look)
      if (col %in% names(df)) m[cbind(si, li)] <- as.numeric(df[[col]])
      m
    }

    enrol_mat <- fill_mat("n.enrolled")
    ev_mat    <- fill_mat("n.event")
    drop_mat  <- fill_mat("n.dropout")
    pipe_mat  <- fill_mat("n.pipeline")
    cut_mat   <- fill_mat("cutoff")

    eff_mat <- matrix(NA_real_, nsim, n_look)
    eff_mat[cbind(si, li)] <- as.numeric(df[[if (use_z) eff.col else p.col]])

    eff_b_mat <- matrix(if (use_z) efficacy else alpha,
                        nsim, n_look, byrow = TRUE)

    # Efficacy crossing (NA statistic or NA boundary -> no crossing)
    reject_cross <- if (use_z && direction == "upper") {
      eff_mat >= eff_b_mat
    } else {
      eff_mat <= eff_b_mat
    }
    reject_cross[is.na(reject_cross)] <- FALSE

    # Futility crossing (Z mode only)
    if (use_fut) {
      fut_mat <- matrix(NA_real_, nsim, n_look)
      fut_mat[cbind(si, li)] <- as.numeric(df[[fut.col]])
      fut_b_mat <- matrix(futility, nsim, n_look, byrow = TRUE)
      fut_cross <- if (direction == "upper") {
        fut_mat <= fut_b_mat
      } else {
        fut_mat >= fut_b_mat
      }
      fut_cross[is.na(fut_cross)] <- FALSE
    } else {
      fut_cross <- matrix(FALSE, nsim, n_look)
    }

    stop_any <- reject_cross | fut_cross

    first_stop <- apply(stop_any, 1L, function(r) {
      w <- which(r)
      if (length(w) > 0L) w[1L] else NA_integer_
    })

    reject   <- logical(nsim)
    reason   <- character(nsim)
    stop_idx <- integer(nsim)
    for (s in seq_len(nsim)) {
      fs <- first_stop[s]
      if (is.na(fs)) {
        reject[s]   <- FALSE
        reason[s]   <- "continue"
        stop_idx[s] <- n_look
      } else {
        stop_idx[s] <- fs
        if (reject_cross[s, fs]) {
          reject[s] <- TRUE
          reason[s] <- "efficacy"
        } else {
          reject[s] <- FALSE
          reason[s] <- "futility"
        }
      }
    }

    rows_idx       <- seq_len(nsim)
    enrol_at_stop  <- enrol_mat[cbind(rows_idx, stop_idx)]
    events_at_stop <- ev_mat[cbind(rows_idx, stop_idx)]
    drop_at_stop   <- drop_mat[cbind(rows_idx, stop_idx)]
    pipe_at_stop   <- pipe_mat[cbind(rows_idx, stop_idx)]
    cutoff_at_stop <- cut_mat[cbind(rows_idx, stop_idx)]

    mean_fin <- function(x) {
      if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
    }

    eff_stop_k <- numeric(n_look)
    fut_stop_k <- numeric(n_look)
    enrol_k    <- numeric(n_look)
    event_k    <- numeric(n_look)
    drop_k     <- numeric(n_look)
    pipe_k     <- numeric(n_look)
    cutoff_k   <- numeric(n_look)
    for (k in seq_len(n_look)) {
      eff_stop_k[k] <- mean(reason == "efficacy" & stop_idx == k)
      fut_stop_k[k] <- mean(reason == "futility" & stop_idx == k)
      enrol_k[k]    <- mean_fin(enrol_mat[, k])
      event_k[k]    <- mean_fin(ev_mat[, k])
      drop_k[k]     <- mean_fin(drop_mat[, k])
      pipe_k[k]     <- mean_fin(pipe_mat[, k])
      cutoff_k[k]   <- mean_fin(cut_mat[, k])
    }

    look_rows <- data.frame(
      population = pop_label,
      look       = as.character(looks),
      stringsAsFactors = FALSE
    )
    if (has_lookval) {
      lv <- tapply(df$look.value, df$look, function(x) x[1L])
      look_rows$look.value <- as.numeric(lv[as.character(looks)])
    }
    if (has_enrol)    look_rows$n.enrolled.mean <- enrol_k
    if (has_event)    look_rows$n.event.mean    <- event_k
    if (has_dropout)  look_rows$n.dropout.mean  <- drop_k
    if (has_pipeline) look_rows$n.pipeline.mean <- pipe_k
    if (has_cutoff)   look_rows$cutoff.mean     <- cutoff_k
    look_rows$prob.stop.efficacy <- eff_stop_k
    look_rows$prob.stop.futility <- fut_stop_k
    look_rows$prob.stop.any      <- eff_stop_k + fut_stop_k
    look_rows$cum.reject         <- cumsum(eff_stop_k)

    overall_row <- data.frame(
      population = pop_label,
      look       = "overall",
      stringsAsFactors = FALSE
    )
    if (has_lookval) overall_row$look.value <- NA_real_
    if (has_enrol)    overall_row$n.enrolled.mean <- mean_fin(enrol_at_stop)
    if (has_event)    overall_row$n.event.mean    <- mean_fin(events_at_stop)
    if (has_dropout)  overall_row$n.dropout.mean  <- mean_fin(drop_at_stop)
    if (has_pipeline) overall_row$n.pipeline.mean <- mean_fin(pipe_at_stop)
    if (has_cutoff)   overall_row$cutoff.mean     <- mean_fin(cutoff_at_stop)
    overall_row$prob.stop.efficacy <- mean(reject)
    overall_row$prob.stop.futility <- mean(reason == "futility")
    overall_row$prob.stop.any      <- mean(reason != "continue")
    overall_row$cum.reject         <- mean(reject)

    list(rows = rbind(look_rows, overall_row), nsim = nsim)
  }

  # ------------------------------------------------------------------ #
  #  Loop over populations and stack
  # ------------------------------------------------------------------ #
  blocks <- vector("list", length(pops))
  nsim_v <- integer(length(pops))
  for (i in seq_along(pops)) {
    df_i <- if (has_pop) {
      data[data$population == pops[i], , drop = FALSE]
    } else {
      data
    }
    r <- summarize_one(df_i, pops[i])
    blocks[[i]] <- r$rows
    nsim_v[i]   <- r$nsim
  }

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL

  structure(out,
            nsim = if (length(unique(nsim_v)) == 1L) nsim_v[1L] else nsim_v,
            boundary = bdesc,
            class = c("simsummary_fast", "data.frame"))
}
