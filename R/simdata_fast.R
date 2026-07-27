#' Fast Simulation of Two-Group Time-to-Event Trial Data
#'
#' @description
#' Simulates time-to-event trial data for one or two groups across many
#' simulated trials, with piecewise accrual, piecewise-exponential survival and
#' dropout, and optional subgroups defined by a prevalence specification. The
#' entire generation pipeline (accrual, survival, dropout, derived columns, and
#' two-group interleaving) runs in a single C++ kernel that materializes the
#' output data frame once, avoiding intermediate R-level vector operations and
#' copies. The random-number stream is consumed in the same order as a per-group
#' reference implementation, so results are reproducible from \code{seed}.
#'
#' @details
#' For each subject the observed time-to-event is
#' \code{tte = pmin(surv_time, dropout_time)} and \code{event} is 1 when the
#' survival time occurs first. The calendar time of the observed event is
#' \code{accrual_time + tte}.
#'
#' The total enrolled is fixed at \code{sum(n)}. With \code{a.rate} the rates are
#' absolute (subjects per unit time): when the accrual period is fully specified
#' the rates must accrue exactly \code{sum(n)}, and when one extra rate is given
#' the end of the final interval is solved so the total is met. With \code{a.prop}
#' the values are relative proportions that distribute \code{sum(n)} across the
#' fully specified intervals. Each accrual interval receives a deterministic
#' number of subjects (the rate or proportion times the group total, rounded to
#' keep the per-group total exact), placed uniformly within the interval.
#'
#' Survival and dropout are exponential when a single hazard (or median) is
#' supplied and piecewise-exponential when a vector is supplied together with
#' the corresponding \code{e.time} or \code{d.time} breakpoints, whose last
#' element must be \code{Inf}. Group-specific parameters are supplied as a
#' two-element list (control first, treatment second).
#'
#' When \code{prevalence} is supplied the trial has subgroups. A numeric vector
#' defines a single factor; a list of numeric vectors defines several
#' independent factors; a multi-dimensional array defines the joint
#' distribution of correlated factors. Per-cell hazards may be supplied as a
#' list with one element per cell. With \code{fixed.alloc = TRUE} the subgroup
#' sizes are deterministic; otherwise subgroup membership is drawn from the
#' prevalence distribution.
#'
#' When \code{n} is a vector of length greater than two together with a per-arm
#' survival list, the simulation is a multi-arm trial. Each arm is generated in
#' turn with the validated single-group kernel over a common accrual window, and
#' the arms are stacked into one data frame with a \code{group} column labeled 1
#' to \code{length(n)} in the order of \code{n}. Per-arm survival is supplied as
#' an \code{e.hazard} or \code{e.median} list with one element per arm, and
#' optional dropout as a shared value or a per-arm list through \code{d.hazard}
#' or \code{d.median}. The arms share the master \code{seed}, so the result is
#' reproducible. A multi-arm design is analyzed as a set of pairwise contrasts by
#' subsetting the output to the control arm and one other arm and calling
#' \code{\link{analysis_fast}} once per contrast. Multi-arm mode does not support
#' subgroups or the illness-death model, which remain two-group.
#'
#' @param nsim Number of simulated trials.
#' @param n Either a single total sample size (split by \code{alloc}), a
#'   length-two vector of per-group sample sizes, or, for a multi-arm trial, a
#'   vector of length greater than two giving the per-arm sample sizes (which
#'   requires a per-arm \code{e.hazard} or \code{e.median} list). When \code{n}
#'   is a per-arm vector, \code{alloc} is ignored.
#' @param alloc A length-two allocation ratio, used when \code{n} is scalar.
#' @param a.time A numeric vector of accrual-interval breakpoints.
#' @param a.rate Absolute accrual rates (subjects per unit time), interpreted in
#'   one of two ways. With length \code{length(a.time) - 1} the accrual period is
#'   fully specified and the rates must accrue exactly \code{sum(n)} subjects (an
#'   inconsistent total is an error). With length \code{length(a.time)} the final
#'   rate applies to an open last interval whose end time is computed so the
#'   total is \code{sum(n)}. Supply exactly one of \code{a.rate} and \code{a.prop}.
#' @param a.prop Accrual proportions, one per accrual interval (length
#'   \code{length(a.time) - 1}), giving the fraction of subjects enrolled in each
#'   interval. Values are normalized to sum to one and distribute the fixed total
#'   \code{sum(n)}. Unlike \code{a.rate} this carries no rate, so the accrual
#'   period must be fully specified by \code{a.time}. Supply exactly one of
#'   \code{a.rate} and \code{a.prop}.
#' @param e.hazard Survival hazard(s). A scalar or vector for one group, or a
#'   two-element list for two groups; per-cell lists are used with subgroups.
#' @param e.median Survival median(s); an alternative to \code{e.hazard}.
#' @param e.time Survival breakpoints for piecewise hazards (last element
#'   \code{Inf}).
#' @param d.hazard Dropout hazard(s), same structure as \code{e.hazard}.
#' @param d.median Dropout median(s); an alternative to \code{d.hazard}.
#' @param d.time Dropout breakpoints for piecewise hazards.
#' @param seed Optional integer seed for the \code{dqrng} generator.
#' @param prevalence Optional subgroup prevalence specification (numeric
#'   vector, list of vectors, array, or a named \code{control}/\code{treatment}
#'   list for group-specific prevalence).
#' @param fixed.alloc Logical; when \code{TRUE} subgroup sizes are
#'   deterministic rather than drawn.
#' @param h01.hazard Transition hazard(s) for the non-terminal (intermediate)
#'   event (state 0 to state 1) in the illness-death model. A scalar or vector
#'   for one group, or a two-element list for two groups. Supplying any of
#'   \code{h01.*} or \code{h02.*} activates the illness-death model with two
#'   correlated endpoints, and is mutually exclusive with \code{e.hazard} /
#'   \code{e.median}.
#' @param h01.median Median(s) for the intermediate event; an alternative to
#'   \code{h01.hazard}.
#' @param h01.time Breakpoints for a piecewise \code{h01.hazard} (last element
#'   \code{Inf}).
#' @param h02.hazard Transition hazard(s) for the terminal event without a prior
#'   intermediate event (state 0 to state 2). Same group and piecewise
#'   conventions as \code{h01.hazard}.
#' @param h02.median Median(s) for the direct terminal event; an alternative to
#'   \code{h02.hazard}.
#' @param h02.time Breakpoints for a piecewise \code{h02.hazard}.
#' @param h12.hazard Transition hazard(s) for the terminal event after an
#'   intermediate event (state 1 to state 2) for subjects who do not switch.
#'   Defaults to \code{h02.hazard}, which gives the Fleischer
#'   maximal-independence model (Fleischer Theorem 1 when there is no switching).
#' @param h12.median Median(s) for the post-event terminal event; an alternative
#'   to \code{h12.hazard}.
#' @param h12.time Breakpoints for a piecewise \code{h12.hazard}, measured from
#'   the intermediate-event time (clock-reset).
#' @param switch.prop Probability that a subject with an intermediate event
#'   switches treatment, a scalar or a two-element list (control, treatment).
#'   Defaults to zero (no switching); the treatment group is typically left at
#'   zero.
#' @param h12.switch.hazard Transition hazard(s) from state 1 to state 2 for
#'   subjects who switch. Required when any \code{switch.prop} is positive.
#' @param h12.switch.median Median(s) for the post-switch terminal event; an
#'   alternative to \code{h12.switch.hazard}.
#' @param h12.switch.time Breakpoints for a piecewise \code{h12.switch.hazard},
#'   measured from the switch (intermediate-event) time (clock-reset).
#' @param switch.clock Time origin for the post-event hazards. Currently only
#'   \code{"reset"} (measured from the intermediate event) is implemented.
#'
#' @return A \code{data.frame} with \code{nsim * sum(n)} rows. The columns are
#'   \code{sim}, \code{group}, any subgroup columns, \code{accrual_time},
#'   \code{surv_time}, \code{dropout_time}, \code{tte}, \code{event}, and
#'   \code{calendar_time}.
#'   In the illness-death model the columns are instead \code{sim},
#'   \code{group}, \code{accrual_time}, \code{e1_surv_time},
#'   \code{e2_surv_time}, \code{dropout_time}, \code{e1_tte}, \code{e1_event},
#'   \code{e2_tte}, \code{e2_event}, \code{e1_calendar_time},
#'   \code{e2_calendar_time}, \code{intermediate}, \code{switched}, and
#'   \code{switch_time}, where \code{e1} is the first (state-0 exit) endpoint
#'   and \code{e2} is the terminal endpoint. In oncology \code{e1} is
#'   progression-free survival, \code{e2} is overall survival, and
#'   \code{intermediate} flags progression.
#'
#' @examples
#' # One-group simulation, simple exponential, no dropout
#' df1 <- simdata_fast(
#'   nsim     = 100,
#'   n        = 50,
#'   a.time   = c(0, 12),
#'   a.rate   = 50 / 12,
#'   e.median = 18,
#'   seed     = 1
#' )
#' head(df1)
#'
#' # Accrual rate with the final interval computed from the total: 20 per unit
#' # time for the first 12 units, then 30 per unit time until 500 are enrolled
#' df1b <- simdata_fast(
#'   nsim     = 100,
#'   n        = 500,
#'   a.time   = c(0, 12),
#'   a.rate   = c(20, 30),
#'   e.median = 18,
#'   seed     = 1
#' )
#' head(df1b)
#'
#' # Accrual by proportion: 30% enrolled in [0, 6], 70% in [6, 12]
#' df1c <- simdata_fast(
#'   nsim     = 100,
#'   n        = 50,
#'   a.time   = c(0, 6, 12),
#'   a.prop   = c(0.3, 0.7),
#'   e.median = 18,
#'   seed     = 1
#' )
#' head(df1c)
#'
#' # Two-group simulation, simple exponential, with dropout
#' df2 <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(60, 60),
#'   a.time   = c(0, 6, 12),
#'   a.rate   = c(8, 12),
#'   e.median = list(18, 24),
#'   d.hazard = list(0.01, 0.01),
#'   seed     = 2
#' )
#' head(df2)
#'
#' # One factor with three levels: single subgroup column
#' df3 <- simdata_fast(
#'   nsim       = 100,
#'   n          = c(150, 150),
#'   a.time     = c(0, 12),
#'   a.rate     = 300 / 12,
#'   e.hazard   = list(list(0.10, 0.08, 0.06), 0.05),
#'   prevalence = c(0.5, 0.3, 0.2),
#'   seed       = 3
#' )
#' head(df3)
#'
#' # Two independent factors (2 x 2): columns subgroup1 and subgroup2.
#' # Four cells in column-major order: (1,1), (2,1), (1,2), (2,2).
#' df4 <- simdata_fast(
#'   nsim       = 100,
#'   n          = 200,
#'   a.time     = c(0, 12),
#'   a.rate     = 200 / 12,
#'   e.hazard   = list(0.10, 0.08, 0.07, 0.05),
#'   prevalence = list(c(0.5, 0.5), c(0.6, 0.4)),
#'   seed       = 4
#' )
#' head(df4)
#'
#' # Two correlated factors via a joint-distribution array (2 x 2)
#' df5 <- simdata_fast(
#'   nsim       = 100,
#'   n          = 200,
#'   a.time     = c(0, 12),
#'   a.rate     = 200 / 12,
#'   e.hazard   = 0.08,
#'   prevalence = array(c(0.40, 0.10, 0.15, 0.35), dim = c(2, 2)),
#'   seed       = 5
#' )
#' head(df5)
#'
#' # Two correlated endpoints, no switching (reduces to Fleischer Theorem 1
#' # because h12 defaults to h02). In oncology e1 is PFS and e2 is OS; here the
#' # control has faster intermediate events and faster direct terminal events.
#' dfid <- simdata_fast(
#'   nsim       = 100,
#'   n          = c(150, 150),
#'   a.time     = c(0, 12),
#'   a.rate     = 300 / 12,
#'   h01.median = list(8, 12),
#'   h02.median = list(24, 30),
#'   seed       = 6
#' )
#' head(dfid)
#'
#' # With treatment switching: 40 percent of control subjects who reach the
#' # intermediate event switch and then follow a more favorable post-event hazard
#' dfsw <- simdata_fast(
#'   nsim              = 100,
#'   n                 = c(150, 150),
#'   a.time            = c(0, 12),
#'   a.rate            = 300 / 12,
#'   h01.median        = list(8, 12),
#'   h02.median        = list(24, 30),
#'   switch.prop       = list(0.4, 0),
#'   h12.switch.median = list(36, 36),
#'   seed              = 7
#' )
#' head(dfsw)
#'
#' # Three-arm trial (one control and two treatment arms) analyzed as pairwise
#' # contrasts against the shared control.
#' dfk <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(120, 120, 120),
#'   a.time   = c(0, 12),
#'   a.rate   = 360 / 12,
#'   e.median = list(12, 16, 20),
#'   seed     = 8
#' )
#' # Control arm (group 1) versus treatment arm 2, one-sided log-rank at month 24.
#' sub12 <- dfk[dfk$group %in% c(1, 2), ]
#' res12 <- analysis_fast(sub12, control = 1, time.looks = 24, side = 1)
#' head(res12)
#'
#' @seealso \code{\link{analysis_fast}}
#'
#' @export
simdata_fast <- function(nsim       = 1000,
                         n,
                         alloc      = c(1, 1),
                         a.time,
                         a.rate     = NULL,
                         a.prop     = NULL,
                         e.hazard   = NULL,
                         e.median   = NULL,
                         e.time     = NULL,
                         d.hazard   = NULL,
                         d.median   = NULL,
                         d.time     = NULL,
                         seed       = NULL,
                         prevalence = NULL,
                         fixed.alloc = FALSE,
                         h01.hazard   = NULL,
                         h01.median   = NULL,
                         h01.time     = NULL,
                         h02.hazard   = NULL,
                         h02.median   = NULL,
                         h02.time     = NULL,
                         h12.hazard   = NULL,
                         h12.median   = NULL,
                         h12.time     = NULL,
                         switch.prop  = NULL,
                         h12.switch.hazard = NULL,
                         h12.switch.median = NULL,
                         h12.switch.time   = NULL,
                         switch.clock = "reset") {

  if (!is.null(seed)) dqrng::dqset.seed(seed)

  # Illness-death (two correlated endpoints, optional switching) mode dispatches
  # to a separate kernel; the single-endpoint path below is left unchanged (same
  # output schema and same dqrng consumption order when these arguments are not
  # supplied).
  if (!is.null(h01.hazard) || !is.null(h01.median) ||
      !is.null(h02.hazard) || !is.null(h02.median)) {
    if (!is.null(e.hazard) || !is.null(e.median)) {
      stop("Specify the illness-death model with 'h01.*' / 'h02.*', not ",
           "'e.hazard' / 'e.median'.")
    }
    return(simdata_fast_id(
      nsim = nsim, n = n, alloc = alloc, a.time = a.time,
      a.rate = a.rate, a.prop = a.prop,
      h01.hazard = h01.hazard, h01.median = h01.median, h01.time = h01.time,
      h02.hazard = h02.hazard, h02.median = h02.median, h02.time = h02.time,
      h12.hazard = h12.hazard, h12.median = h12.median, h12.time = h12.time,
      switch.prop = switch.prop,
      h12.switch.hazard = h12.switch.hazard,
      h12.switch.median = h12.switch.median,
      h12.switch.time = h12.switch.time,
      switch.clock = switch.clock,
      d.hazard = d.hazard, d.median = d.median, d.time = d.time,
      prevalence = prevalence))
  }

  # Multi-arm (K > 2) mode dispatches to a separate assembler that generates each
  # arm with the validated single-group kernel and stacks the results with
  # per-arm group labels 1 to K. It is entered only when 'n' has length greater
  # than two AND a matching per-arm survival list is supplied; a length-(>2) 'n'
  # without such a list falls through to the historical length check below, so
  # the one-group and two-group paths, their dqrng consumption order, and the
  # existing input-validation behavior are all unchanged. The master seed set
  # above is shared by every arm, so the multi-arm result is reproducible.
  if (length(n) > 2L) {
    surv_spec <- if (!is.null(e.hazard)) e.hazard else e.median
    if (is.list(surv_spec) && length(surv_spec) == length(n)) {
      return(simdata_fast_karm(
        nsim = nsim, n = n, a.time = a.time,
        a.rate = a.rate, a.prop = a.prop,
        e.hazard = e.hazard, e.median = e.median, e.time = e.time,
        d.hazard = d.hazard, d.median = d.median, d.time = d.time,
        prevalence = prevalence,
        h01.hazard = h01.hazard, h01.median = h01.median,
        h02.hazard = h02.hazard, h02.median = h02.median))
    }
  }

  use_subgroup        <- !is.null(prevalence)
  group_specific_prev <- use_subgroup && is_group_specific_prev(prevalence)

  spec_ctrl <- NULL
  spec_trt  <- NULL
  if (use_subgroup) {
    if (group_specific_prev) {
      spec_ctrl <- build_prevalence_spec(prevalence[["control"]])
      spec_trt  <- build_prevalence_spec(prevalence[["treatment"]])
      if (spec_ctrl$n_cell != spec_trt$n_cell ||
          !identical(dim(spec_ctrl$level_table), dim(spec_trt$level_table))) {
        stop("'control' and 'treatment' prevalence must use the same ",
             "number of factors and levels")
      }
    } else {
      spec_ctrl <- build_prevalence_spec(prevalence)
      spec_trt  <- spec_ctrl
    }
  }

  if (length(n) == 1L) {
    n_grp    <- round(n * alloc / sum(alloc))
    n_groups <- 2L
  } else if (length(n) == 2L) {
    n_grp    <- n
    n_groups <- 2L
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }

  if (use_subgroup) {
    is_two_group <- (length(n) == 2L) || group_specific_prev
  } else {
    is_two_group <- n_groups == 2L && (
      is.list(e.hazard) || is.list(e.median) ||
        is.list(d.hazard) || is.list(d.median)
    )
    if (!is_two_group && length(n) == 2L) {
      stop("Two-group sample sizes supplied but hazard/median parameters are not lists.\n",
           "Wrap group-specific parameters in list().")
    }
  }

  if (!is_two_group) {
    n_grp    <- if (length(n) == 1L) n else n[1L]
    n_groups <- 1L
  }

  if (use_subgroup && is_two_group) {
    chk_two <- function(x, nm) {
      if (!is.null(x) && !(is.list(x) && length(x) == 2L)) {
        stop("For a two-group simulation, '", nm,
             "' must be a list of length 2 (one element per group).")
      }
    }
    chk_two(e.hazard, "e.hazard")
    chk_two(e.median, "e.median")
    chk_two(d.hazard, "d.hazard")
    chk_two(d.median, "d.median")
  }

  if (!is.null(e.hazard) && !is.null(e.median)) {
    stop("Specify exactly one of 'e.hazard' and 'e.median', not both")
  }
  if (is.null(e.hazard) && is.null(e.median)) {
    stop("One of 'e.hazard' or 'e.median' must be supplied")
  }
  if (!is.null(e.median)) e.hazard <- convert_median_to_hazard(e.median)

  has_dropout <- !is.null(d.hazard) || !is.null(d.median)
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify exactly one of 'd.hazard' and 'd.median', not both")
  }
  if (!is.null(d.median)) d.hazard <- convert_median_to_hazard(d.median)

  use_rate <- !is.null(a.rate)
  use_prop <- !is.null(a.prop)
  if (use_rate == use_prop) {
    stop("Supply exactly one of 'a.rate' and 'a.prop'")
  }
  if (length(a.time) < 2L) stop("'a.time' must have at least two elements")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  n_total    <- sum(n_grp)
  n_int_time <- length(a.time) - 1L
  acc_tol    <- 1e-8 * max(1, n_total)

  if (use_rate) {
    if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
    if (length(a.rate) == n_int_time) {
      implied <- sum(a.rate * diff(a.time))
      if (abs(implied - n_total) > acc_tol) {
        stop("'a.rate' implies ", round(implied, 4), " subjects over the ",
             "accrual period but sum(n) is ", n_total, ".\n",
             "Make them consistent, drop the final 'a.time' breakpoint to let ",
             "it be computed from the rate, or use 'a.prop' for relative ",
             "accrual.")
      }
      a.time_full <- a.time
    } else if (length(a.rate) == n_int_time + 1L) {
      bounded <- if (n_int_time >= 1L) {
        sum(a.rate[seq_len(n_int_time)] * diff(a.time))
      } else {
        0
      }
      remaining <- n_total - bounded
      if (remaining <= 0) {
        stop("The specified accrual intervals already accrue ",
             round(bounded, 4), " subjects, at least sum(n) = ", n_total,
             "; the final interval cannot be extended.")
      }
      final_dur   <- remaining / a.rate[length(a.rate)]
      a.time_full <- c(a.time, a.time[length(a.time)] + final_dur)
    } else {
      stop("'a.rate' must have length equal to length(a.time) - 1 (fully ",
           "specified accrual) or length(a.time) (final interval end computed ",
           "from sum(n))")
    }
    weights_a <- a.rate * diff(a.time_full)
  } else {
    if (any(a.prop <= 0)) stop("All 'a.prop' values must be positive")
    if (length(a.prop) != n_int_time) {
      stop("'a.prop' must have length equal to length(a.time) - 1")
    }
    a.time_full <- a.time
    weights_a   <- as.numeric(a.prop)
  }

  a_int_prob <- weights_a / sum(weights_a)

  n_cell_ctrl <- if (use_subgroup) spec_ctrl$n_cell else 1L
  n_cell_trt  <- if (use_subgroup) spec_trt$n_cell  else 1L
  n_cell      <- n_cell_ctrl

  # Build the per-cell hazard / piecewise pre-computation lists for one group.
  # Each returned list has one element per cell: the hazard vector, the finite
  # breakpoints, and the cumulative hazard at those breakpoints. For a single
  # hazard the breakpoint / cumulative entries are empty (unused by the kernel).
  build_exp_specs <- function(hazard_grp, time_grp, n_cell) {
    haz <- vector("list", n_cell)
    fin <- vector("list", n_cell)
    cum <- vector("list", n_cell)
    for (c in seq_len(n_cell)) {
      hz <- if (is.list(hazard_grp)) hazard_grp[[c]] else hazard_grp
      tm <- if (n_cell > 1L && is.list(time_grp)) time_grp[[c]] else time_grp
      pc <- piecewise_precompute(hz, tm)
      haz[[c]] <- pc$hazard
      fin[[c]] <- pc$fin_time
      cum[[c]] <- pc$cum_haz
    }
    list(haz = haz, fin = fin, cum = cum)
  }

  # Resolve group-level survival hazard specs for control and treatment.
  if (use_subgroup && is_two_group) {
    e_haz_c <- if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard
    e_haz_t <- if (is.list(e.hazard)) e.hazard[[2L]] else e.hazard
  } else if (n_groups == 2L && !use_subgroup) {
    e_haz_c <- if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard
    e_haz_t <- if (is.list(e.hazard)) e.hazard[[2L]] else e.hazard
  } else {
    e_haz_c <- e.hazard
    e_haz_t <- e.hazard
  }

  d_haz_c <- NULL; d_haz_t <- NULL
  if (has_dropout) {
    if (use_subgroup && is_two_group) {
      d_haz_c <- if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard
      d_haz_t <- if (is.list(d.hazard)) d.hazard[[2L]] else d.hazard
    } else if (n_groups == 2L && !use_subgroup) {
      d_haz_c <- if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard
      d_haz_t <- if (is.list(d.hazard)) d.hazard[[2L]] else d.hazard
    } else {
      d_haz_c <- if (use_subgroup) d.hazard else if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard
      d_haz_t <- d_haz_c
    }
  }

  e.time_c <- resolve_time_arg(e.time, 1L)
  e.time_t <- resolve_time_arg(e.time, 2L)
  d.time_c <- if (has_dropout) resolve_time_arg(d.time, 1L) else NULL
  d.time_t <- if (has_dropout) resolve_time_arg(d.time, 2L) else NULL

  e_c <- build_exp_specs(e_haz_c, e.time_c, n_cell_ctrl)
  e_t <- build_exp_specs(e_haz_t, e.time_t, n_cell_trt)
  d_c <- if (has_dropout) build_exp_specs(d_haz_c, d.time_c, n_cell_ctrl) else
    list(haz = list(numeric(0)), fin = list(numeric(0)), cum = list(numeric(0)))
  d_t <- if (has_dropout) build_exp_specs(d_haz_t, d.time_t, n_cell_trt) else
    list(haz = list(numeric(0)), fin = list(numeric(0)), cum = list(numeric(0)))

  # Subgroup descriptors.
  if (use_subgroup) {
    cum_prev_c   <- spec_ctrl$cum_prev
    cum_prev_t   <- spec_trt$cum_prev
    level_tab_c  <- matrix(as.integer(spec_ctrl$level_table),
                           nrow = nrow(spec_ctrl$level_table))
    level_tab_t  <- matrix(as.integer(spec_trt$level_table),
                           nrow = nrow(spec_trt$level_table))
    sub_names    <- spec_ctrl$sub_names
    fixed_c      <- if (fixed.alloc) fixed_cell_counts(n_grp[1L], spec_ctrl$cell_prob) else integer(n_cell_ctrl)
    fixed_t      <- if (fixed.alloc) fixed_cell_counts(if (n_groups == 2L) n_grp[2L] else n_grp[1L], spec_trt$cell_prob) else integer(n_cell_trt)
  } else {
    cum_prev_c <- numeric(0); cum_prev_t <- numeric(0)
    level_tab_c <- matrix(integer(0), nrow = 1L, ncol = 0L)
    level_tab_t <- matrix(integer(0), nrow = 1L, ncol = 0L)
    sub_names  <- character(0)
    fixed_c    <- integer(0); fixed_t <- integer(0)
  }

  n_grp_int <- if (n_groups == 1L) as.integer(n_grp[1L]) else as.integer(n_grp[1:2])

  # Deterministic per-interval accrual counts for each group. Each group's
  # counts sum to its size exactly (largest-remainder rounding), and the kernel
  # repeats them every simulation, placing subjects uniformly within intervals.
  acc_counts_c <- accrual_cell_counts(n_grp_int[1L], a_int_prob)
  acc_counts_t <- if (n_groups == 2L) {
    accrual_cell_counts(n_grp_int[2L], a_int_prob)
  } else {
    integer(0)
  }

  simdata_core_full(
    as.integer(nsim), n_grp_int,
    as.numeric(a.time_full), acc_counts_c, acc_counts_t,
    as.integer(n_cell),
    e_c$haz, e_c$fin, e_c$cum,
    e_t$haz, e_t$fin, e_t$cum,
    has_dropout,
    d_c$haz, d_c$fin, d_c$cum,
    d_t$haz, d_t$fin, d_t$cum,
    as.numeric(cum_prev_c), as.numeric(cum_prev_t),
    level_tab_c, level_tab_t,
    as.character(sub_names),
    fixed.alloc,
    as.integer(fixed_c), as.integer(fixed_t)
  )
}

# Piecewise-exponential pre-computation: validates the hazard / breakpoint
# pair and returns the hazard vector, the finite left breakpoints, and the
# cumulative hazard at those breakpoints, matching rpiece_exp_r. For a single
# hazard the breakpoint / cumulative entries are empty.
piecewise_precompute <- function(hazard, e.time) {
  if (length(hazard) == 1L) {
    return(list(hazard = as.numeric(hazard),
                fin_time = numeric(0), cum_haz = numeric(0)))
  }
  if (is.null(e.time)) {
    stop("'e.time' must be supplied when 'e.hazard' (or 'e.median') is a vector")
  }
  if (length(e.time) != length(hazard) + 1L) {
    stop("length(e.time) must equal length(e.hazard) + 1")
  }
  if (!is.infinite(e.time[length(e.time)])) {
    stop("Last element of 'e.time' must be Inf")
  }
  n_int    <- length(hazard)
  fin_time <- e.time[seq_len(n_int)]
  lengths  <- diff(e.time[seq_len(n_int)])
  cum_haz  <- c(0, cumsum(hazard[-n_int] * lengths))
  list(hazard = as.numeric(hazard),
       fin_time = as.numeric(fin_time), cum_haz = as.numeric(cum_haz))
}

# ------------------------------------------------------------------ #
#  Internal helpers: subgroup prevalence and parameter resolution
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
#  Internal helper: detect group-specific prevalence
# ------------------------------------------------------------------ #
is_group_specific_prev <- function(prev) {
  is.list(prev) && !is.array(prev) && !is.null(names(prev)) &&
    setequal(names(prev), c("control", "treatment"))
}

# ------------------------------------------------------------------ #
#  Internal helper: build a prevalence specification for one group
# ------------------------------------------------------------------ #
# Returns a list with the cumulative cell probabilities (cum_prev), a
# cell-by-factor level table (level_table, n_cell x n_fac), the output column
# names (sub_names), the number of cells (n_cell), and the number of factors
# (n_fac). Cells are ordered column-major (first factor varies fastest).
build_prevalence_spec <- function(prev) {
  if (is.array(prev) && length(dim(prev)) >= 2L) {
    # Joint distribution of correlated factors
    dims <- dim(prev)
    cell_prob <- as.numeric(prev)
    if (any(cell_prob <= 0)) stop("All 'prevalence' values must be positive")
    cell_prob   <- cell_prob / sum(cell_prob)
    level_table <- as.matrix(expand.grid(lapply(dims, seq_len)))
    n_fac       <- length(dims)
  } else if (is.list(prev)) {
    # Several independent factors (marginal prevalence per factor)
    bad <- vapply(prev, function(x) !is.numeric(x) || is.list(x), logical(1L))
    if (any(bad)) {
      stop("Each factor in a 'prevalence' list must be a numeric vector")
    }
    if (any(unlist(prev) <= 0)) stop("All 'prevalence' values must be positive")
    margins     <- lapply(prev, function(x) x / sum(x))
    levels_f    <- vapply(margins, length, integer(1L))
    level_table <- as.matrix(expand.grid(lapply(levels_f, seq_len)))
    cell_prob   <- apply(level_table, 1L, function(idx) {
      prod(mapply(function(m, i) m[i], margins, idx))
    })
    n_fac <- length(margins)
  } else {
    # Single factor
    if (!is.numeric(prev)) stop("'prevalence' must be numeric")
    if (any(prev <= 0)) stop("All 'prevalence' values must be positive")
    cell_prob   <- prev / sum(prev)
    level_table <- matrix(seq_along(prev), ncol = 1L)
    n_fac       <- 1L
  }

  n_cell    <- length(cell_prob)
  sub_names <- if (n_fac == 1L) "subgroup" else paste0("subgroup", seq_len(n_fac))
  dimnames(level_table) <- NULL

  list(cum_prev = cumsum(cell_prob), cell_prob = cell_prob,
       level_table = level_table, sub_names = sub_names,
       n_cell = n_cell, n_fac = n_fac)
}

# ------------------------------------------------------------------ #
#  Internal helper: deterministic cell counts for fixed allocation
# ------------------------------------------------------------------ #
# Each cell receives floor(n * p) patients; the remaining patients (n minus
# the sum of the floors) are added one at a time starting from the first
# cell, so the last cells absorb the rounding shortfall.
fixed_cell_counts <- function(n, p) {
  base <- floor(n * p)
  rem  <- n - sum(base)
  if (rem > 0L) base[seq_len(rem)] <- base[seq_len(rem)] + 1L
  as.integer(base)
}

# ------------------------------------------------------------------ #
#  Internal helper: largest-remainder accrual counts per interval
# ------------------------------------------------------------------ #
# Distributes n subjects across accrual intervals in proportion to p, with the
# rounding remainder assigned to the intervals with the largest fractional
# parts (Hamilton's method). This keeps the per-interval counts as close as
# possible to n * p and is robust to floating-point error in p, so that e.g.
# proportions of c(0.48, 0.52) with n = 500 give exactly c(240, 260).
accrual_cell_counts <- function(n, p) {
  target <- n * p
  base   <- floor(target)
  rem    <- n - sum(base)
  if (rem > 0L) {
    take <- order(target - base, decreasing = TRUE)[seq_len(rem)]
    base[take] <- base[take] + 1L
  }
  as.integer(base)
}

# ------------------------------------------------------------------ #
#  Internal helper: resolve a possibly group-specific time argument
# ------------------------------------------------------------------ #
resolve_time_arg <- function(time_arg, group_idx) {
  if (is.null(time_arg)) return(NULL)
  if (is.list(time_arg)) time_arg[[group_idx]] else time_arg
}

# ------------------------------------------------------------------ #
#  Internal helper: convert median survival time to hazard
# ------------------------------------------------------------------ #
convert_median_to_hazard <- function(median_arg) {
  if (is.list(median_arg)) {
    lapply(median_arg, convert_median_to_hazard)
  } else {
    log(2) / median_arg
  }
}


# ------------------------------------------------------------------ #
#  Internal helper: resolve accrual interval counts (illness-death path)
# ------------------------------------------------------------------ #
# Replicates the accrual resolution of the single-endpoint path for the
# illness-death wrapper, returning the completed accrual breakpoints and the
# per-interval subject counts for each group.
resolve_accrual_counts <- function(n_grp_int, a.time, a.rate, a.prop) {
  n_total    <- sum(n_grp_int)
  n_int_time <- length(a.time) - 1L
  acc_tol    <- 1e-8 * max(1, n_total)
  use_rate   <- !is.null(a.rate)
  use_prop   <- !is.null(a.prop)
  if (use_rate == use_prop) stop("Supply exactly one of 'a.rate' and 'a.prop'")
  if (length(a.time) < 2L) stop("'a.time' must have at least two elements")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  if (use_rate) {
    if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
    if (length(a.rate) == n_int_time) {
      implied <- sum(a.rate * diff(a.time))
      if (abs(implied - n_total) > acc_tol) {
        stop("'a.rate' implies ", round(implied, 4), " subjects over the ",
             "accrual period but sum(n) is ", n_total, ".")
      }
      a.time_full <- a.time
    } else if (length(a.rate) == n_int_time + 1L) {
      bounded <- if (n_int_time >= 1L) {
        sum(a.rate[seq_len(n_int_time)] * diff(a.time))
      } else {
        0
      }
      remaining <- n_total - bounded
      if (remaining <= 0) {
        stop("The specified accrual intervals already accrue ", round(bounded, 4),
             " subjects, at least sum(n) = ", n_total, ".")
      }
      final_dur   <- remaining / a.rate[length(a.rate)]
      a.time_full <- c(a.time, a.time[length(a.time)] + final_dur)
    } else {
      stop("'a.rate' must have length equal to length(a.time) - 1 or length(a.time)")
    }
    weights_a <- a.rate * diff(a.time_full)
  } else {
    if (any(a.prop <= 0)) stop("All 'a.prop' values must be positive")
    if (length(a.prop) != n_int_time) {
      stop("'a.prop' must have length equal to length(a.time) - 1")
    }
    a.time_full <- a.time
    weights_a   <- as.numeric(a.prop)
  }

  a_int_prob   <- weights_a / sum(weights_a)
  acc_counts_c <- accrual_cell_counts(as.integer(n_grp_int[1L]), a_int_prob)
  acc_counts_t <- if (length(n_grp_int) == 2L) {
    accrual_cell_counts(as.integer(n_grp_int[2L]), a_int_prob)
  } else {
    integer(0)
  }
  list(a.time_full = a.time_full,
       acc_counts_c = acc_counts_c, acc_counts_t = acc_counts_t)
}

# ------------------------------------------------------------------ #
#  Internal wrapper: illness-death (two-endpoint, switching) simulation
# ------------------------------------------------------------------ #
# Validates the transition-hazard arguments, builds the per-group
# piecewise-exponential pre-computations, and calls simdata_core_id. The seed is
# already set by the caller. Subgroups are not supported here.
simdata_fast_id <- function(nsim, n, alloc, a.time, a.rate, a.prop,
                            h01.hazard, h01.median, h01.time,
                            h02.hazard, h02.median, h02.time,
                            h12.hazard, h12.median, h12.time,
                            switch.prop,
                            h12.switch.hazard, h12.switch.median, h12.switch.time,
                            switch.clock,
                            d.hazard, d.median, d.time,
                            prevalence) {

  if (!is.null(prevalence)) {
    stop("Subgroups ('prevalence') are not supported with the illness-death ",
         "model; call simdata_fast separately per subgroup and combine with rbind().")
  }
  if (!identical(switch.clock, "reset")) {
    stop("Only switch.clock = \"reset\" is currently implemented.")
  }

  # Resolve transition hazards (exactly one of hazard or median per quantity).
  if (!is.null(h01.hazard) && !is.null(h01.median)) {
    stop("Specify exactly one of 'h01.hazard' and 'h01.median'")
  }
  if (is.null(h01.hazard) && is.null(h01.median)) {
    stop("One of 'h01.hazard' or 'h01.median' must be supplied for the ",
         "illness-death model")
  }
  if (!is.null(h01.median)) h01.hazard <- convert_median_to_hazard(h01.median)

  if (!is.null(h02.hazard) && !is.null(h02.median)) {
    stop("Specify exactly one of 'h02.hazard' and 'h02.median'")
  }
  if (is.null(h02.hazard) && is.null(h02.median)) {
    stop("One of 'h02.hazard' or 'h02.median' must be supplied for the ",
         "illness-death model")
  }
  if (!is.null(h02.median)) h02.hazard <- convert_median_to_hazard(h02.median)

  if (!is.null(h12.hazard) && !is.null(h12.median)) {
    stop("Specify exactly one of 'h12.hazard' and 'h12.median'")
  }
  if (!is.null(h12.median)) h12.hazard <- convert_median_to_hazard(h12.median)

  if (!is.null(h12.switch.hazard) && !is.null(h12.switch.median)) {
    stop("Specify exactly one of 'h12.switch.hazard' and 'h12.switch.median'")
  }
  if (!is.null(h12.switch.median)) {
    h12.switch.hazard <- convert_median_to_hazard(h12.switch.median)
  }

  has_dropout <- !is.null(d.hazard) || !is.null(d.median)
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify exactly one of 'd.hazard' and 'd.median'")
  }
  if (!is.null(d.median)) d.hazard <- convert_median_to_hazard(d.median)

  # A two-group simulation is signalled by a length-two 'n' or any group-specific
  # (length-two list) hazard / switch argument.
  is_two_group <- (length(n) == 2L) ||
    is.list(h01.hazard) || is.list(h02.hazard) ||
    is.list(h12.hazard) || is.list(h12.switch.hazard) ||
    is.list(switch.prop)

  if (length(n) == 1L) {
    n_grp <- if (is_two_group) round(n * alloc / sum(alloc)) else n
  } else if (length(n) == 2L) {
    n_grp <- n
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }
  n_groups  <- if (is_two_group) 2L else 1L
  n_grp_int <- if (n_groups == 1L) as.integer(n_grp[1L]) else as.integer(n_grp[1:2])
  g2        <- (n_groups == 2L)

  # Resolve a possibly group-specific argument for group g (1 = control,
  # 2 = treatment). A length-two list is per group; anything else is shared.
  id_grp <- function(x, g) if (is.list(x)) x[[g]] else x

  get_pi <- function(g) {
    if (is.null(switch.prop)) return(0)
    p <- id_grp(switch.prop, g)
    if (is.null(p)) 0 else as.numeric(p)
  }
  pi_c <- get_pi(1L)
  pi_t <- if (g2) get_pi(2L) else 0

  if ((pi_c > 0 || pi_t > 0) && is.null(h12.switch.hazard)) {
    stop("'h12.switch.hazard' (or 'h12.switch.median') must be supplied when ",
         "any 'switch.prop' is positive")
  }

  # Post-event no-switch hazard defaults to the direct terminal hazard h02.
  get_h12      <- function(g) if (is.null(h12.hazard)) id_grp(h02.hazard, g) else id_grp(h12.hazard, g)
  get_h12_time <- function(g) if (is.null(h12.hazard)) id_grp(h02.time, g)   else id_grp(h12.time, g)

  h01_c <- piecewise_precompute(id_grp(h01.hazard, 1L), id_grp(h01.time, 1L))
  h01_t <- if (g2) piecewise_precompute(id_grp(h01.hazard, 2L), id_grp(h01.time, 2L)) else h01_c
  h02_c <- piecewise_precompute(id_grp(h02.hazard, 1L), id_grp(h02.time, 1L))
  h02_t <- if (g2) piecewise_precompute(id_grp(h02.hazard, 2L), id_grp(h02.time, 2L)) else h02_c
  h12_c <- piecewise_precompute(get_h12(1L), get_h12_time(1L))
  h12_t <- if (g2) piecewise_precompute(get_h12(2L), get_h12_time(2L)) else h12_c

  empty_spec <- list(hazard = 1, fin_time = numeric(0), cum_haz = numeric(0))
  if (pi_c > 0 || pi_t > 0) {
    h12s_c <- piecewise_precompute(id_grp(h12.switch.hazard, 1L), id_grp(h12.switch.time, 1L))
    h12s_t <- if (g2) {
      piecewise_precompute(id_grp(h12.switch.hazard, 2L), id_grp(h12.switch.time, 2L))
    } else {
      h12s_c
    }
  } else {
    h12s_c <- empty_spec
    h12s_t <- empty_spec
  }

  if (has_dropout) {
    d_c <- piecewise_precompute(id_grp(d.hazard, 1L), id_grp(d.time, 1L))
    d_t <- if (g2) piecewise_precompute(id_grp(d.hazard, 2L), id_grp(d.time, 2L)) else d_c
  } else {
    d_c <- empty_spec
    d_t <- empty_spec
  }

  acc <- resolve_accrual_counts(n_grp_int, a.time, a.rate, a.prop)

  simdata_core_id(
    as.integer(nsim), n_grp_int,
    as.numeric(acc$a.time_full), acc$acc_counts_c, acc$acc_counts_t,
    h01_c$hazard, h01_c$fin_time, h01_c$cum_haz,
    h01_t$hazard, h01_t$fin_time, h01_t$cum_haz,
    h02_c$hazard, h02_c$fin_time, h02_c$cum_haz,
    h02_t$hazard, h02_t$fin_time, h02_t$cum_haz,
    h12_c$hazard, h12_c$fin_time, h12_c$cum_haz,
    h12_t$hazard, h12_t$fin_time, h12_t$cum_haz,
    h12s_c$hazard, h12s_c$fin_time, h12s_c$cum_haz,
    h12s_t$hazard, h12s_t$fin_time, h12s_t$cum_haz,
    as.numeric(pi_c), as.numeric(pi_t),
    has_dropout,
    d_c$hazard, d_c$fin_time, d_c$cum_haz,
    d_t$hazard, d_t$fin_time, d_t$cum_haz
  )
}

# ------------------------------------------------------------------ #
#  Internal helper: resolve whole-trial accrual for the multi-arm path
# ------------------------------------------------------------------ #
# Completes the accrual breakpoints (solving an open final interval from the
# total when a trailing 'a.rate' is supplied) and returns the per-interval
# accrual proportions, so every arm can be generated over a common accrual
# window by passing the completed breakpoints together with these proportions as
# a fully specified 'a.prop'. This mirrors the accrual resolution of the
# single-endpoint path and is deliberately self-contained so the two-group path
# is left byte-identical.
resolve_karm_accrual <- function(n_total, a.time, a.rate, a.prop) {
  n_int_time <- length(a.time) - 1L
  acc_tol    <- 1e-8 * max(1, n_total)
  use_rate   <- !is.null(a.rate)
  use_prop   <- !is.null(a.prop)
  if (use_rate == use_prop) stop("Supply exactly one of 'a.rate' and 'a.prop'")
  if (length(a.time) < 2L) stop("'a.time' must have at least two elements")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  if (use_rate) {
    if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
    if (length(a.rate) == n_int_time) {
      implied <- sum(a.rate * diff(a.time))
      if (abs(implied - n_total) > acc_tol) {
        stop("'a.rate' implies ", round(implied, 4), " subjects over the ",
             "accrual period but sum(n) is ", n_total, ".")
      }
      a.time_full <- a.time
    } else if (length(a.rate) == n_int_time + 1L) {
      bounded <- if (n_int_time >= 1L) {
        sum(a.rate[seq_len(n_int_time)] * diff(a.time))
      } else {
        0
      }
      remaining <- n_total - bounded
      if (remaining <= 0) {
        stop("The specified accrual intervals already accrue ", round(bounded, 4),
             " subjects, at least sum(n) = ", n_total, ".")
      }
      final_dur   <- remaining / a.rate[length(a.rate)]
      a.time_full <- c(a.time, a.time[length(a.time)] + final_dur)
    } else {
      stop("'a.rate' must have length equal to length(a.time) - 1 or length(a.time)")
    }
    weights_a <- a.rate * diff(a.time_full)
  } else {
    if (any(a.prop <= 0)) stop("All 'a.prop' values must be positive")
    if (length(a.prop) != n_int_time) {
      stop("'a.prop' must have length equal to length(a.time) - 1")
    }
    a.time_full <- a.time
    weights_a   <- as.numeric(a.prop)
  }
  list(a.time_full = a.time_full, a_int_prob = weights_a / sum(weights_a))
}

# ------------------------------------------------------------------ #
#  Internal wrapper: multi-arm (K > 2) simulation
# ------------------------------------------------------------------ #
# Generates each of the K arms with the validated single-group path of
# simdata_fast (a scalar 'n' and a non-list hazard), sharing a common accrual
# window supplied as fully specified proportions, and stacks the arms into one
# data frame with a 'group' column labeled 1 to K in the order of 'n'. The master
# seed is set by the caller and the arms consume the same dqrng stream in order,
# so the result is reproducible from 'seed'. Subgroups and the illness-death
# model are not supported here.
simdata_fast_karm <- function(nsim, n, a.time, a.rate, a.prop,
                              e.hazard, e.median, e.time,
                              d.hazard, d.median, d.time,
                              prevalence,
                              h01.hazard, h01.median,
                              h02.hazard, h02.median) {
  K <- length(n)

  if (!is.null(prevalence)) {
    stop("Multi-arm mode (length(n) > 2) does not support 'prevalence' ",
         "subgroups.")
  }
  if (!is.null(h01.hazard) || !is.null(h01.median) ||
      !is.null(h02.hazard) || !is.null(h02.median)) {
    stop("Multi-arm mode (length(n) > 2) does not support the illness-death ",
         "model, which is two-group.")
  }
  if (!is.null(e.hazard) && !is.null(e.median)) {
    stop("Specify exactly one of 'e.hazard' and 'e.median'")
  }
  if (is.null(e.hazard) && is.null(e.median)) {
    stop("One of 'e.hazard' or 'e.median' must be supplied")
  }

  # Per-arm survival specification: a length-K list, one element per arm.
  surv_is_median <- is.null(e.hazard)
  surv_arg       <- if (surv_is_median) e.median else e.hazard
  if (!is.list(surv_arg) || length(surv_arg) != K) {
    stop("In multi-arm mode '", if (surv_is_median) "e.median" else "e.hazard",
         "' must be a list with one element per arm, so its length must equal ",
         "length(n) = ", K, ".")
  }

  # Optional per-arm dropout: NULL, a shared scalar or vector, or a length-K list.
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify at most one of 'd.hazard' and 'd.median'")
  }
  drop_is_median <- is.null(d.hazard) && !is.null(d.median)
  drop_arg       <- if (!is.null(d.hazard)) d.hazard else d.median
  has_dropout    <- !is.null(drop_arg)
  if (has_dropout && is.list(drop_arg) && length(drop_arg) != K) {
    stop("In multi-arm mode a list 'd.hazard' or 'd.median' must have one ",
         "element per arm, so its length must equal length(n) = ", K, ".")
  }

  # Resolve the whole-trial accrual once so every arm shares the same accrual
  # window and shape.
  acc         <- resolve_karm_accrual(sum(n), a.time, a.rate, a.prop)
  a.time_full <- acc$a.time_full
  a_int_prob  <- acc$a_int_prob

  # A possibly per-arm 'e.time' / 'd.time': a list is indexed per arm, anything
  # else is shared across arms.
  arm_e_time <- function(g) if (is.list(e.time)) e.time[[g]] else e.time
  arm_d_time <- function(g) if (is.list(d.time)) d.time[[g]] else d.time

  blocks <- vector("list", K)
  for (g in seq_len(K)) {
    args <- list(
      nsim   = nsim,
      n      = as.numeric(n[g]),
      a.time = a.time_full,
      a.prop = a_int_prob,
      e.time = arm_e_time(g)
    )
    if (surv_is_median) {
      args$e.median <- surv_arg[[g]]
    } else {
      args$e.hazard <- surv_arg[[g]]
    }
    if (has_dropout) {
      dv <- if (is.list(drop_arg)) drop_arg[[g]] else drop_arg
      if (drop_is_median) args$d.median <- dv else args$d.hazard <- dv
      args$d.time <- arm_d_time(g)
    }
    dfg         <- do.call(simdata_fast, args)
    dfg$group   <- g
    blocks[[g]] <- dfg
  }

  out <- do.call(rbind, blocks)
  # Group rows by simulation (then by arm) so the output matches the
  # single-endpoint contract that the analysis path relies on for its skip-sort
  # fast path.
  out <- out[order(out$sim, out$group), , drop = FALSE]
  rownames(out) <- NULL
  out
}
