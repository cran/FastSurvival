#' Fast Treatment Switching in Simulated Trial Data
#'
#' @description
#' Applies a treatment-switching rule to simulated trial data, such as the
#' output of \code{\link{simdata_fast}}, by changing only the outcomes that
#' occur after each subject's switch. Subjects of the selected groups switch
#' either at their intermediate event (for example disease progression), at
#' an opening time after an interim analysis (crossover at a milestone), or at
#' the later of the two, and the remaining time to the terminal event after the
#' switch is either multiplied by an acceleration factor or redrawn from a new
#' (piecewise) exponential hazard. All simulated trials are processed at once
#' with vectorized operations. Because the outcomes observed before the switch
#' are kept exactly, an analysis at or before the opening time gives the same
#' result before and after switching, so a crossover that depends on an interim
#' result can be simulated by analyzing the interim, deciding per simulated
#' trial, and switching afterward.
#'
#' @details
#' The function works on the latent columns of \code{\link{simdata_fast}}. For
#' a single-endpoint simulation these are \code{surv_time} and
#' \code{dropout_time}; for an illness-death simulation they are
#' \code{e1_surv_time} (time to the first event, such as progression-free
#' survival), \code{e2_surv_time} (time to the terminal event, such as overall
#' survival), \code{dropout_time}, and \code{intermediate}. The observed
#' columns (\code{tte}, \code{event}, and \code{calendar_time}, or their
#' \code{e1_} and \code{e2_} counterparts) are recomputed from the modified
#' latent times, and the columns \code{switched} and \code{switch_time} (time
#' from accrual to the switch) record the switches. The columns \code{sim},
#' \code{group}, \code{accrual_time}, and the subgroup columns are unchanged.
#'
#' The switch time \code{s}, measured from accrual, is
#' \itemize{
#'   \item \code{when = "intermediate"}: the intermediate event time
#'     \code{e1_surv_time} of a subject whose intermediate event occurred
#'     (\code{intermediate == 1}); illness-death data only;
#'   \item \code{when = "cutoff"}: the opening time \code{cutoff + delay} of the
#'     subject's simulated trial minus the accrual time (zero for a subject
#'     accrued after the opening);
#'   \item \code{when = "later"}: the later of the intermediate event time and
#'     the opening time, for a subject whose intermediate event occurred;
#'     illness-death data only.
#' }
#' A subject of the selected \code{group} switches with probability
#' \code{prob} when the switch occurs before the terminal event and before
#' dropout (\code{s < e2_surv_time}, or \code{s < surv_time} for
#' single-endpoint data, and \code{s < dropout_time}) and the subject has not
#' already switched (\code{switched == 1}, for example at progression through
#' the \code{switch.prop} argument of \code{\link{simdata_fast}}). With \code{when = "cutoff"} or \code{"later"}, the switch
#' never occurs before the opening time, so every outcome observed by calendar
#' time \code{cutoff + delay} is left unchanged; subjects accrued after the
#' opening are also eligible, as in a protocol that opens crossover from that
#' time on.
#'
#' The effect of switching on the terminal time \code{T > s} is
#' \itemize{
#'   \item \code{aft.factor = f}: \code{T_new = s + f * (T - s)}, the causal
#'     accelerated failure time model of the rank-preserving structural failure
#'     time (RPSFT) method, under which \code{f > 1} prolongs the remaining
#'     survival;
#'   \item \code{hazard} (or \code{median}): \code{T_new = s + R}, where the
#'     remaining time \code{R} is drawn from the (piecewise) exponential
#'     distribution with hazard \code{hazard} and breakpoints \code{time}
#'     measured from the switch.
#' }
#' In illness-death data with \code{when = "cutoff"}, a switch can occur before
#' the intermediate event. The acceleration factor is then applied in the same
#' way to the first-event time when it lies after the switch, which keeps the
#' first event no later than the terminal event; a redrawn hazard is not
#' available in this case because the latent time to the intermediate event
#' after the terminal event is not stored.
#'
#' \code{sims} restricts the switching to some simulated trials, for example
#' those in which the interim analysis crossed an efficacy boundary. It is
#' available with \code{when = "cutoff"} or \code{"later"} only, so that the
#' decision cannot change outcomes observed before it is made.
#'
#' The random numbers (one uniform per row when \code{prob < 1}, and one
#' exponential per row when \code{hazard} or \code{median} is used) are drawn
#' with \code{dqrng}, in row order, so the result is reproducible from
#' \code{seed} and \code{stream}, or from the state left by a seeded
#' \code{\link{simdata_fast}} call.
#'
#' @param data A data frame from \code{\link{simdata_fast}} (single-endpoint,
#'   multi-arm, or illness-death), containing at least \code{sim},
#'   \code{group}, \code{accrual_time}, and the latent columns described in
#'   Details.
#' @param group A vector of group labels whose subjects may switch, for example
#'   the control group \code{1}.
#' @param prob A single probability in [0, 1] that an eligible subject
#'   switches. Defaults to 1.
#' @param when A character string, one of \code{"intermediate"},
#'   \code{"cutoff"}, and \code{"later"}, giving the switch time (see Details).
#' @param cutoff Per-simulation calendar times at which switching opens, used
#'   with \code{when = "cutoff"} or \code{"later"}. Either a numeric vector with
#'   one element per simulated trial or a one-column matrix from
#'   \code{\link{cutoff_fast}}. As in the \code{cutoff.looks} argument of
#'   \code{\link{analysis_fast}}, the names of a vector (or the row names of a
#'   matrix) are matched to the values of \code{data$sim}; without names the
#'   elements are taken in the order of the sorted distinct values of
#'   \code{data$sim}. \code{NA} disables switching in that simulated trial.
#' @param delay A single non-negative number added to \code{cutoff}, for
#'   example the time needed to implement the crossover after an interim
#'   analysis. Defaults to 0.
#' @param sims An optional logical vector with one element per simulated trial
#'   (in the order of the sorted distinct values of \code{data$sim}). Switching
#'   is applied only in the simulated trials with \code{TRUE}. Requires
#'   \code{when = "cutoff"} or \code{"later"}.
#' @param aft.factor A single positive acceleration factor applied to the
#'   remaining time after the switch. Exactly one of \code{aft.factor},
#'   \code{hazard}, and \code{median} must be supplied.
#' @param hazard Post-switch hazard, a single value or a vector for a piecewise
#'   exponential distribution with breakpoints \code{time}, measured from the
#'   switch.
#' @param median A single post-switch median; an alternative to a single
#'   \code{hazard}.
#' @param time Breakpoints for a piecewise \code{hazard}, starting at 0 and
#'   ending with \code{Inf}, measured from the switch.
#' @param seed Optional integer seed for the \code{dqrng} generator. If
#'   \code{NULL} (default), the draws come from the current state of that
#'   generator, which \code{set.seed()} does not control; supply \code{seed}
#'   for reproducible results (see \code{\link{simdata_fast}}).
#' @param stream Optional non-negative whole number selecting a \code{dqrng}
#'   stream; requires \code{seed}. See \code{\link{simdata_fast}}.
#'
#' @return The data frame \code{data} with the latent and observed columns of
#'   the switched subjects updated, and the columns \code{switched} (0 or 1)
#'   and \code{switch_time} (time from accrual to the switch, \code{NA} for
#'   subjects who did not switch) added or updated.
#'
#' @references
#' Robins, J. M., & Tsiatis, A. A. (1991). Correcting for non-compliance in
#' randomized trials using rank preserving structural failure time models.
#' \emph{Communications in Statistics - Theory and Methods}, \emph{20}(8),
#' 2609-2631.
#'
#' @examples
#' # Crossover at an interim analysis: control subjects still on study at the
#' # 150th event switch to the experimental treatment, which multiplies their
#' # remaining survival time by 1.5.
#' df <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
#'                    a.rate = 300 / 12, e.median = list(12, 18),
#'                    d.hazard = 0.01, seed = 1)
#' ia <- cutoff_fast(df, event.looks = 150)
#' dfs <- switch_fast(df, group = 1, when = "cutoff", cutoff = ia,
#'                    aft.factor = 1.5)
#' mean(dfs$switched[dfs$group == 1])
#'
#' # The interim analysis is unchanged; the final analysis is diluted.
#' r0 <- analysis_fast(df,  control = 1, cutoff.looks = ia)
#' r1 <- analysis_fast(dfs, control = 1, cutoff.looks = ia)
#' all.equal(r0$logrank.z, r1$logrank.z)
#' fa0 <- analysis_fast(df,  control = 1, event.looks = 220)
#' fa1 <- analysis_fast(dfs, control = 1, event.looks = 220)
#' c(mean(fa0$logrank.z), mean(fa1$logrank.z))
#'
#' # Switching at progression in an illness-death simulation: 50 percent of
#' # the control subjects who progress switch, with a post-progression median
#' # of 15 instead of 10.
#' id <- simdata_fast(nsim = 50, n = c(150, 150), a.time = c(0, 12),
#'                    a.rate = 300 / 12, h01.median = list(6, 9),
#'                    h02.median = list(20, 24), h12.median = list(10, 10),
#'                    seed = 2)
#' ids <- switch_fast(id, group = 1, prob = 0.5, when = "intermediate",
#'                    median = 15)
#' table(ids$group, ids$switched)
#'
#' @seealso \code{\link{simdata_fast}}, \code{\link{cutoff_fast}},
#'   \code{\link{analysis_fast}}.
#'
#' @export
switch_fast <- function(data, group, prob = 1,
                        when = c("intermediate", "cutoff", "later"),
                        cutoff = NULL, delay = 0, sims = NULL,
                        aft.factor = NULL, hazard = NULL, median = NULL,
                        time = NULL, seed = NULL, stream = NULL) {
  when <- match.arg(when)

  # ---- Data type ---------------------------------------------------------
  if (!is.data.frame(data) ||
      !all(c("sim", "group", "accrual_time") %in% names(data))) {
    stop("'data' must be a data frame with columns 'sim', 'group', and ",
         "'accrual_time', as produced by simdata_fast()")
  }
  id_cols <- c("e1_surv_time", "e2_surv_time", "dropout_time", "intermediate")
  is_id <- all(id_cols %in% names(data))
  if (!is_id && !all(c("surv_time", "dropout_time") %in% names(data))) {
    stop("'data' must contain the latent columns of simdata_fast(): ",
         "'surv_time' and 'dropout_time', or the illness-death columns ",
         paste(id_cols, collapse = ", "))
  }
  if (!is_id && when != "cutoff") {
    stop("when = \"", when, "\" needs an intermediate event; use ",
         "when = \"cutoff\" for single-endpoint data")
  }

  # ---- Arguments ---------------------------------------------------------
  if (missing(group) || length(group) < 1L || !all(group %in% data$group)) {
    stop("'group' must give one or more group labels present in 'data'")
  }
  if (length(prob) != 1L || !is.finite(prob) || prob < 0 || prob > 1) {
    stop("'prob' must be a single value in [0, 1]")
  }
  if (length(delay) != 1L || !is.finite(delay) || delay < 0) {
    stop("'delay' must be a single non-negative value")
  }
  n_eff <- (!is.null(aft.factor)) + (!is.null(hazard)) + (!is.null(median))
  if (n_eff != 1L) {
    stop("supply exactly one of 'aft.factor', 'hazard', and 'median'")
  }
  use_aft <- !is.null(aft.factor)
  if (use_aft) {
    if (length(aft.factor) != 1L || !is.finite(aft.factor) ||
        aft.factor <= 0) {
      stop("'aft.factor' must be a single positive value")
    }
    if (!is.null(time)) stop("'time' is used only with 'hazard'")
  } else {
    if (!is.null(median)) {
      if (length(median) != 1L || !is.finite(median) || median <= 0) {
        stop("'median' must be a single positive value")
      }
      if (!is.null(time)) stop("'time' is used only with 'hazard'")
      hazard <- log(2) / median
    }
    if (!is.numeric(hazard) || length(hazard) < 1L || anyNA(hazard) ||
        any(hazard < 0) || any(is.infinite(hazard)) ||
        hazard[length(hazard)] <= 0) {
      stop("'hazard' must be finite and non-negative, with a positive last ",
           "value")
    }
    if (length(hazard) == 1L && !is.null(time)) {
      stop("'time' is used only with a piecewise (vector) 'hazard'")
    }
    if (length(hazard) > 1L) {
      if (is.null(time) || length(time) != length(hazard) + 1L ||
          anyNA(time) || time[1L] != 0 || any(diff(time) <= 0) ||
          !is.infinite(time[length(time)])) {
        stop("'time' must start at 0, be strictly increasing, end with Inf, ",
             "and have one more element than 'hazard'")
      }
    }
    if (is_id && when == "cutoff") {
      stop("with illness-death data and when = \"cutoff\", only 'aft.factor' ",
           "is available (see Details)")
    }
    pc <- piecewise_precompute(hazard, time)
  }

  sim_ids <- sort(unique(data$sim))
  nsim    <- length(sim_ids)
  row_sim <- match(data$sim, sim_ids)

  uses_open <- when %in% c("cutoff", "later")
  if (uses_open) {
    if (is.null(cutoff)) {
      stop("'cutoff' must be supplied when when = \"", when, "\"")
    }
    if (!is.matrix(cutoff) && !is.null(names(cutoff))) {
      cutoff <- matrix(as.numeric(cutoff), ncol = 1L,
                       dimnames = list(names(cutoff), NULL))
    }
    if (is.matrix(cutoff)) {
      if (ncol(cutoff) != 1L) {
        stop("a matrix 'cutoff' must have one column; select the look, for ",
             "example cutoff[, 1, drop = FALSE]")
      }
      rn <- rownames(cutoff)
      cv <- as.numeric(cutoff[, 1L])
      if (!is.null(rn)) {
        idx <- match(sim_key(sim_ids), rn)
        if (anyNA(idx)) {
          stop("the names (or row names) of 'cutoff' do not cover every ",
               "value of 'data$sim'")
        }
        cv <- cv[idx]
      }
    } else {
      cv <- as.numeric(cutoff)
    }
    if (length(cv) != nsim) {
      stop("'cutoff' must have one element per simulated trial (", nsim, ")")
    }
    if (any(!is.na(cv) & (cv < 0 | is.infinite(cv)))) {
      stop("'cutoff' must be non-negative and finite (or NA)")
    }
    open_sim <- cv + delay
  } else {
    if (!is.null(cutoff)) {
      stop("'cutoff' is used only with when = \"cutoff\" or \"later\"")
    }
    open_sim <- rep(0, nsim)
  }
  if (!is.null(sims)) {
    if (!uses_open) {
      stop("'sims' requires when = \"cutoff\" or \"later\", so that the ",
           "switch follows the decision")
    }
    if (!is.logical(sims) || length(sims) != nsim || anyNA(sims)) {
      stop("'sims' must be a logical vector with one non-missing element per ",
           "simulated trial (", nsim, ")")
    }
  } else {
    sims <- rep(TRUE, nsim)
  }

  if (!is.null(stream) && is.null(seed)) stop("'stream' requires 'seed'")
  if (!is.null(seed)) {
    if (!is.null(stream) &&
        (length(stream) != 1L || !is.finite(stream) || stream < 0 ||
         stream > .Machine$integer.max ||
         abs(stream - round(stream)) > 1e-8)) {
      stop("'stream' must be a single non-negative whole number")
    }
    # dqrng would truncate a value such as 100 * 0.29 to the integer below.
    if (!is.null(stream)) stream <- round(stream)
    dqrng::dqset.seed(seed, stream = stream)
  }

  # ---- Random numbers (row order) -----------------------------------------
  n_row <- nrow(data)
  u <- if (prob < 1) dqrng::dqrunif(n_row) else NULL
  e <- if (!use_aft) dqrng::dqrexp(n_row, rate = 1) else NULL

  # ---- Switch time and eligibility ---------------------------------------
  a      <- as.numeric(data$accrual_time)
  drop_t <- as.numeric(data$dropout_time)
  term   <- if (is_id) as.numeric(data$e2_surv_time) else as.numeric(data$surv_time)
  open_c <- open_sim[row_sim]
  open_r <- pmax(open_c - a, 0)

  s <- switch(when,
              intermediate = as.numeric(data$e1_surv_time),
              cutoff       = open_r,
              later        = pmax(as.numeric(data$e1_surv_time), open_r))
  already <- if ("switched" %in% names(data)) data$switched %in% 1 else
    rep(FALSE, n_row)
  elig <- data$group %in% group & sims[row_sim] & !already &
    !is.na(s) & s < term & s < drop_t
  if (when != "cutoff") elig <- elig & data$intermediate %in% 1
  if (uses_open) {
    # Compare on the calendar scale with the same sum (accrual + time) that
    # cutoff_fast() and analysis_fast() use, so that a subject whose event or
    # dropout is at the opening time (for example the event that defines an
    # event-driven cutoff) is never switched because of rounding in
    # open_c - a.
    elig <- elig & (a + term > open_c) & (a + drop_t > open_c)
  }
  if (prob < 1) elig <- elig & u < prob
  idx <- which(elig)

  # ---- Post-switch terminal time -----------------------------------------
  s_i <- s[idx]
  if (use_aft) {
    term_new <- s_i + aft.factor * (term[idx] - s_i)
  } else {
    term_new <- s_i + inv_piecewise(e[idx], pc)
  }
  if (uses_open) {
    # A modified time must stay after the opening time on the calendar scale;
    # otherwise (possible only through rounding) the original time is kept.
    keep <- !(a[idx] + term_new > open_c[idx])
    term_new[keep] <- term[idx][keep]
  }

  out <- data
  sw_col <- if ("switched" %in% names(out)) as.integer(out$switched) else
    integer(n_row)
  st_col <- if ("switch_time" %in% names(out)) as.numeric(out$switch_time) else
    rep(NA_real_, n_row)
  sw_col[idx] <- 1L
  st_col[idx] <- s_i

  if (is_id) {
    e1 <- as.numeric(out$e1_surv_time)
    e2 <- as.numeric(out$e2_surv_time)
    if (when == "cutoff") {
      # The first event after the switch is accelerated in the same way.
      late <- idx[a[idx] + e1[idx] > open_c[idx]]
      s_l  <- s[late]
      e1_new <- s_l + aft.factor * (e1[late] - s_l)
      ok_l <- a[late] + e1_new > open_c[late]
      e1[late[ok_l]] <- e1_new[ok_l]
    }
    e2[idx] <- term_new
    out$e1_surv_time     <- e1
    out$e2_surv_time     <- e2
    # Ties count as events; a subject with neither a finite event time nor a
    # finite dropout time is censored, as in simdata_fast().
    out$e1_tte           <- pmin(e1, drop_t)
    out$e1_event         <- as.integer(e1 < drop_t |
                                         (e1 == drop_t & is.finite(e1)))
    out$e2_tte           <- pmin(e2, drop_t)
    out$e2_event         <- as.integer(e2 < drop_t |
                                         (e2 == drop_t & is.finite(e2)))
    out$e1_calendar_time <- a + out$e1_tte
    out$e2_calendar_time <- a + out$e2_tte
  } else {
    st <- as.numeric(out$surv_time)
    st[idx] <- term_new
    out$surv_time <- st
    out$tte   <- pmin(st, drop_t)
    out$event <- as.integer(st < drop_t | (st == drop_t & is.finite(st)))
    if ("calendar_time" %in% names(out)) out$calendar_time <- a + out$tte
  }
  out$switched    <- sw_col
  out$switch_time <- st_col
  out
}

# ------------------------------------------------------------------ #
#  Internal helper: inverse cumulative hazard of a piecewise exponential
# ------------------------------------------------------------------ #
# 'pc' is the output of piecewise_precompute(); 'e' are unit exponential
# draws. The time t solves H(t) = e, and is Inf when the last hazard is zero
# and e exceeds the cumulative hazard at the last breakpoint.
inv_piecewise <- function(e, pc) {
  if (length(pc$hazard) == 1L) {
    return(if (pc$hazard > 0) e / pc$hazard else rep(Inf, length(e)))
  }
  k <- findInterval(e, pc$cum_haz)
  h <- pc$hazard[k]
  ifelse(h > 0, pc$fin_time[k] + (e - pc$cum_haz[k]) / h, Inf)
}
