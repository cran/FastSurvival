#' Fast Per-Simulation Analysis Cutoffs from Combined Trigger Rules
#'
#' @description
#' Computes the calendar time of each analysis look in every simulated trial
#' from a combination of trigger rules: a target number of events, a planned
#' calendar time, a maximum calendar time, a minimum time after the previous
#' look, and a minimum follow-up after a given number of subjects are enrolled.
#' The cutoffs are computed for all simulated trials in a single C++ pass and
#' are returned as a matrix that can be passed to the \code{cutoff.looks}
#' argument of \code{\link{analysis_fast}} and \code{\link{pairwise_fast}}, or
#' to the \code{cutoff} argument of \code{\link{switch_fast}}. Separating the
#' determination of the analysis time from the analysis itself lets the same
#' cutoffs drive several endpoints (for example overall survival analyzed at
#' the progression-free survival event-driven looks), several contrasts, or a
#' modification of the data after an interim analysis.
#'
#' @details
#' For look \code{l} of a simulated trial the candidate calendar times are
#' \itemize{
#'   \item \code{t_cal}, the planned calendar time \code{time.looks[l]};
#'   \item \code{t_event}, the calendar time (accrual time plus observed time)
#'     of the \code{d}-th counted event, where \code{d} is \code{event.looks[l]}
#'     (or the entry of the \code{event.looks} matrix for that simulation and
#'     look);
#'   \item \code{t_gap}, the cutoff of the previous look plus \code{min.gap[l]}
#'     (zero plus \code{min.gap[1]} at the first look);
#'   \item \code{t_enr}, the accrual time of the \code{min.enrolled[l]}-th
#'     enrolled subject plus \code{min.followup[l]};
#' }
#' and the cutoff is
#' \code{min(max(t_cal, t_event, t_gap, t_enr), max.time[l])}, where a
#' condition that is not supplied (or is \code{NA} at that look) is left out.
#' This is the rule of \code{get_analysis_date()} in the simtrial package
#' (with \code{max.time} playing the role of its
#' \code{max_extension_for_target_event}), applied to every simulated trial at
#' once. For example, \code{event.looks = 300} with \code{max.time = 48} is "300
#' events or month 48, whichever comes first", and \code{event.looks = 300}
#' with \code{time.looks = 36} is "300 events but not before month 36".
#'
#' An event or enrollment target that is not met in the simulated data has no
#' finite calendar time. Unless a \code{max.time} caps the look, its cutoff is
#' then \code{NA}, which \code{\link{analysis_fast}} reports as
#' \code{reached = FALSE} with the statistics of the full data, as it does for
#' an unreached \code{event.looks} target. In this respect the function differs
#' from simtrial, which uses the last observed event time for an unmet target.
#'
#' Events are counted on the columns named by \code{tte.col} and
#' \code{event.col}, so a look can be triggered by an endpoint other than the
#' one analyzed (for example the \code{e1_tte} and \code{e1_event} columns of
#' an illness-death simulation). By default every subject's events are
#' counted; \code{event.subset} restricts the count to a subset of rows, such as
#' the control group (\code{data$group == 1}), one subgroup, or the two arms of
#' the primary contrast in a multi-arm trial. When the subjects tie at the
#' target event time, all events at that time are included in the analysis,
#' so the analyzed event count can exceed the target, as in
#' \code{\link{analysis_fast}}.
#'
#' A matrix \code{event.looks} with one row per simulated trial gives
#' simulation-specific event targets, for example targets re-estimated at an
#' interim analysis.
#'
#' Each look is determined by its own conditions; only \code{min.gap} refers
#' to the previous look. A look with \code{max.time} as its only condition is
#' placed at that time. When the previous look is not reached, \code{t_gap} is
#' infinite, so a look with \code{min.gap} is not reached either unless its
#' \code{max.time} caps it. The looks are not reordered: when the conditions
#' of a later look give an earlier cutoff than the previous look (for example
#' \code{event.looks = c(150, 220)} with \code{max.time = c(NA, 30)}), a
#' warning is given. Supplying \code{min.gap} (for example 0) at the later
#' looks keeps them in order unless their \code{max.time} is earlier.
#'
#' @param data A data frame of simulated trials, such as the output of
#'   \code{\link{simdata_fast}}, with columns \code{sim}, \code{accrual_time},
#'   and the columns named by \code{tte.col} and \code{event.col}.
#' @param event.looks Target cumulative event counts. Either a vector with one
#'   positive whole number per look (\code{NA} to omit the event condition at a
#'   look), or a matrix of positive whole numbers with one row per simulated
#'   trial (in the order of the sorted distinct values of \code{data$sim}) and
#'   one column per look.
#' @param time.looks Planned calendar times, one per look (\code{NA} to omit).
#' @param max.time Maximum calendar times, one per look (\code{NA} to omit).
#'   The cutoff never exceeds this value.
#' @param min.gap Minimum time after the previous look, one per look
#'   (\code{NA} to omit). At the first look it is measured from time zero.
#' @param min.enrolled Number of enrolled subjects after which the minimum
#'   follow-up \code{min.followup} must elapse, one per look (\code{NA} to
#'   omit).
#' @param min.followup Minimum follow-up after the \code{min.enrolled}-th
#'   subject is enrolled, one per look. Defaults to zero when
#'   \code{min.enrolled} is supplied, and requires \code{min.enrolled}.
#' @param event.subset An optional logical vector with one element per row of
#'   \code{data}. When supplied, only the events of the rows with \code{TRUE}
#'   are counted for \code{event.looks}.
#' @param tte.col A single character string naming the observed-time column
#'   (numeric and non-negative) used to count events. Defaults to
#'   \code{"tte"}.
#' @param event.col A single character string naming the event-indicator column
#'   (numeric or logical, 0 or 1) used to count events. Defaults to
#'   \code{"event"}.
#'
#' @return A numeric matrix with one row per simulated trial and one column per
#'   look, holding the calendar cutoffs (\code{NA} for a look that is not
#'   reached). The row names are the simulation identifiers, and the attribute
#'   \code{"look.value"} holds, for each look, the event target when a common
#'   target is supplied, otherwise the planned calendar time, otherwise
#'   \code{NA}. \code{\link{analysis_fast}} reports this attribute in its
#'   \code{look.value} column.
#'
#' @examples
#' df <- simdata_fast(nsim = 20, n = c(150, 150), a.time = c(0, 12),
#'                    a.rate = 300 / 12, e.median = list(12, 18),
#'                    d.hazard = 0.01, seed = 1)
#'
#' # 150 events or month 30, whichever comes first, then a final look at
#' # 220 events but at least 6 months after the interim
#' cut <- cutoff_fast(df, event.looks = c(150, 220), max.time = c(30, NA),
#'                    min.gap = c(NA, 6))
#' head(cut)
#'
#' res <- analysis_fast(df, control = 1, cutoff.looks = cut)
#' head(res)
#'
#' # Looks triggered by the events of the control group only
#' cut_c <- cutoff_fast(df, event.looks = 100, event.subset = df$group == 1)
#' head(cut_c)
#'
#' @seealso \code{\link{analysis_fast}}, \code{\link{pairwise_fast}},
#'   \code{\link{switch_fast}}, \code{\link{simdata_fast}}.
#'
#' @export
cutoff_fast <- function(data, event.looks = NULL, time.looks = NULL,
                        max.time = NULL, min.gap = NULL,
                        min.enrolled = NULL, min.followup = NULL,
                        event.subset = NULL,
                        tte.col = "tte", event.col = "event") {

  # ---- Columns -----------------------------------------------------------
  if (!is.character(tte.col) || length(tte.col) != 1L ||
      !is.character(event.col) || length(event.col) != 1L) {
    stop("'tte.col' and 'event.col' must be single character strings")
  }
  req_cols <- c("sim", "accrual_time", tte.col, event.col)
  if (!is.data.frame(data) || !all(req_cols %in% names(data))) {
    stop("'data' must be a data frame with columns: ",
         paste(req_cols, collapse = ", "))
  }
  if (nrow(data) == 0L) stop("'data' has no rows")
  if (anyNA(data$sim) || anyNA(data$accrual_time) || anyNA(data[[tte.col]])) {
    stop("columns 'sim', 'accrual_time', and '", tte.col, "' of 'data' must ",
         "not contain missing values")
  }
  if (!is.numeric(data[[tte.col]]) || any(data[[tte.col]] < 0)) {
    stop("column '", tte.col, "' of 'data' must be numeric and non-negative")
  }
  ev <- data[[event.col]]
  if (!(is.numeric(ev) || is.logical(ev))) {
    stop("column '", event.col, "' of 'data' must be numeric or logical, not ",
         "a factor or a character vector")
  }
  if (anyNA(ev) || !all(ev == 0 | ev == 1)) {
    stop("column '", event.col, "' of 'data' must be coded as 0 or 1")
  }
  if (!is.null(event.subset)) {
    if (!is.logical(event.subset) || length(event.subset) != nrow(data) ||
        anyNA(event.subset)) {
      stop("'event.subset' must be a logical vector with one non-missing ",
           "element per row of 'data'")
    }
  }

  # ---- Simulation offsets ------------------------------------------------
  sim_vec <- data$sim
  if (is.unsorted(sim_vec)) {
    ord   <- order(sim_vec)
    reidx <- function(x) x[ord]
  } else {
    reidx <- function(x) x
  }
  sim_s   <- reidx(sim_vec)
  sim_ids <- unique(sim_s)
  nsim    <- length(sim_ids)
  counts  <- tabulate(match(sim_s, sim_ids), nbins = nsim)
  sim_ptr <- as.integer(c(0L, cumsum(counts)))

  # ---- Number of looks ---------------------------------------------------
  ev_is_mat <- is.matrix(event.looks)
  lens <- c(
    if (ev_is_mat) ncol(event.looks) else length(event.looks),
    length(time.looks), length(max.time), length(min.gap),
    length(min.enrolled), length(min.followup)
  )
  n_looks <- max(lens)
  if (n_looks < 1L) {
    stop("supply at least one of 'event.looks', 'time.looks', 'max.time', ",
         "'min.gap', and 'min.enrolled'")
  }
  bad_len <- lens > 1L & lens != n_looks
  if (any(bad_len)) {
    stop("each look specification must have length 1 or the number of looks (",
         n_looks, ")")
  }

  # Expand a per-look specification to length n_looks (NA when absent).
  expand <- function(x) {
    if (is.null(x)) return(rep(NA_real_, n_looks))
    x <- as.numeric(x)
    if (length(x) == 1L) rep(x, n_looks) else x
  }
  t_cal <- expand(time.looks)
  t_cap <- expand(max.time)
  gap   <- expand(min.gap)
  n_enr <- expand(min.enrolled)
  f_up  <- expand(min.followup)

  check_pos <- function(x, nm, allow_zero = FALSE) {
    y <- x[!is.na(x)]
    if (any(!is.finite(y)) || any(if (allow_zero) y < 0 else y <= 0)) {
      stop("'", nm, "' must be ", if (allow_zero) "non-negative" else "positive",
           " and finite (or NA)")
    }
  }
  check_pos(t_cal, "time.looks")
  check_pos(t_cap, "max.time")
  check_pos(gap, "min.gap", allow_zero = TRUE)
  check_pos(f_up, "min.followup", allow_zero = TRUE)
  check_pos(n_enr, "min.enrolled")
  if (any(!is.na(n_enr) & abs(n_enr - round(n_enr)) > 1e-8)) {
    stop("'min.enrolled' must be positive whole numbers (or NA)")
  }
  n_enr <- round(n_enr)
  if (!is.null(min.followup) && is.null(min.enrolled)) {
    stop("'min.followup' requires 'min.enrolled'")
  }
  f_up[!is.na(n_enr) & is.na(f_up)] <- 0

  # Event targets as an nsim x n_looks matrix (NA when absent).
  if (is.null(event.looks)) {
    target <- matrix(NA_real_, nrow = nsim, ncol = n_looks)
  } else if (ev_is_mat) {
    if (nrow(event.looks) != nsim) {
      stop("a matrix 'event.looks' must have one row per simulated trial (",
           nsim, ")")
    }
    tg <- as.numeric(event.looks)
    if (anyNA(tg) || any(!is.finite(tg)) || any(tg < 1) ||
        any(abs(tg - round(tg)) > 1e-8)) {
      stop("a matrix 'event.looks' must contain positive whole numbers")
    }
    target <- matrix(round(tg), nrow = nsim, ncol = n_looks)
  } else {
    tg <- expand(event.looks)
    y  <- tg[!is.na(tg)]
    if (any(!is.finite(y)) || any(y < 1) || any(abs(y - round(y)) > 1e-8)) {
      stop("'event.looks' must be positive whole numbers (or NA)")
    }
    target <- matrix(round(tg), nrow = nsim, ncol = n_looks, byrow = TRUE)
  }

  # Targets beyond the integer range can never be met; cap them so that the
  # C++ conversion to int is defined.
  int_max <- .Machine$integer.max
  target[!is.na(target) & target > int_max] <- int_max
  n_enr[!is.na(n_enr) & n_enr > int_max] <- int_max

  has_lower <- !is.na(t_cal) | !is.na(gap) | !is.na(n_enr) |
    colSums(!is.na(target)) > 0
  if (any(!has_lower & is.na(t_cap))) {
    stop("every look needs at least one trigger condition")
  }

  count_ev <- if (is.null(event.subset)) {
    rep(1L, nrow(data))
  } else {
    as.integer(event.subset)
  }

  core <- cutoff_core(
    sim_ptr,
    reidx(as.numeric(data$accrual_time)),
    reidx(as.numeric(data[[tte.col]])),
    reidx(as.integer(ev)),
    reidx(count_ev),
    t_cal, target, t_cap, gap, n_enr, f_up
  )
  core[!is.finite(core)] <- NA_real_

  # Each look follows its own rules, so a later look can fall before the
  # previous one (for example when only the later look has a calendar cap).
  # The looks are not reordered; the user is warned instead.
  if (n_looks > 1L) {
    back   <- core[, -1L, drop = FALSE] < core[, -n_looks, drop = FALSE]
    n_back <- sum(rowSums(back, na.rm = TRUE) > 0)
    if (n_back > 0L) {
      warning("in ", n_back, " simulated trial(s) the cutoff of a look is ",
              "earlier than that of the previous look; see Details of ",
              "?cutoff_fast", call. = FALSE)
    }
  }

  look_value <- if (!is.null(event.looks) && !ev_is_mat) {
    ifelse(!is.na(target[1L, ]), target[1L, ], t_cal)
  } else {
    t_cal
  }
  dimnames(core) <- list(sim_key(sim_ids), paste0("look", seq_len(n_looks)))
  attr(core, "look.value") <- as.numeric(look_value)
  core
}
