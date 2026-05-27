#' Fast Survival Data Simulation for Clinical Trials
#'
#' @description
#' Simulates individual patient data for one- or two-group time-to-event
#' trials. Patient accrual times are drawn from a piecewise uniform
#' distribution. Survival and dropout times are drawn from either a simple
#' exponential or a piecewise exponential distribution, depending on whether
#' a scalar or vector hazard is supplied. Random number generation uses
#' \code{dqrng} for speed. C++ backends handle piecewise sampling and
#' two-group interleaving to minimize R-level overhead.
#'
#' @details
#' For each patient, three times are generated independently:
#'
#' \describe{
#'   \item{Accrual time}{Drawn from a piecewise uniform distribution defined
#'     by \code{a.time} and \code{a.rate}. The probability of falling in
#'     interval j is proportional to \code{a.rate[j] * (a.time[j+1] - a.time[j])}.}
#'   \item{Survival time}{Drawn from a simple exponential distribution when
#'     \code{e.hazard} is a scalar (or \code{e.median} is a scalar), or from
#'     a piecewise exponential distribution when \code{e.hazard} is a vector
#'     (with \code{e.time} required). The inverse-CDF method is used for the
#'     piecewise case.}
#'   \item{Dropout time}{Drawn from the same family as survival time, governed
#'     by \code{d.hazard} / \code{d.median} and \code{d.time}. Set to
#'     \code{Inf} when \code{d.hazard} and \code{d.median} are both
#'     \code{NULL} (no dropout).}
#' }
#'
#' The observed time-to-event is \code{tte = pmin(surv_time, dropout_time)}.
#' The event indicator is \code{1} when \code{surv_time <= dropout_time} and
#' \code{0} otherwise. The calendar time is
#' \code{calendar_time = accrual_time + tte}.
#'
#' Exactly one of \code{e.hazard} and \code{e.median} must be supplied.
#' Exactly one of \code{d.hazard} and \code{d.median} must be supplied when
#' dropout is modeled; both must be \code{NULL} to suppress dropout entirely.
#'
#' For a two-group trial, group-specific parameters are passed as a list of
#' length 2, where element 1 corresponds to the control group and element 2
#' to the treatment group. For a one-group trial, plain vectors or scalars
#' are passed directly.
#'
#' @param nsim A positive integer specifying the number of simulation
#'   iterations. Default is 1000.
#' @param n A positive integer (one group) or a numeric vector of length 2
#'   (two groups: \code{c(n_control, n_treatment)}) specifying the per-group
#'   sample sizes per simulation. Alternatively, supply the total sample size
#'   as a scalar together with \code{alloc} to split automatically.
#' @param alloc A numeric vector of length 2 specifying the allocation ratio
#'   \code{c(control, treatment)}. Used only when \code{n} is a scalar.
#'   Default is \code{c(1, 1)}.
#' @param a.time A numeric vector of accrual interval boundaries
#'   \code{c(0, t1, ..., T)} with \code{a.time[1] = 0}. Must be strictly
#'   increasing. The last element is the end of the accrual window.
#' @param a.rate A numeric vector of accrual rates (patients per time unit)
#'   for each interval. Length must equal \code{length(a.time) - 1}. All
#'   values must be positive.
#' @param e.hazard Hazard rate(s) for the event time. For a simple exponential,
#'   supply a scalar (one group) or a list of two scalars (two groups). For a
#'   piecewise exponential, supply a numeric vector (one group) or a list of
#'   two numeric vectors (two groups); \code{e.time} must then also be
#'   supplied. Cannot be used together with \code{e.median}.
#' @param e.median Median survival time(s) for the event time. Same structure
#'   as \code{e.hazard}. Converted internally via
#'   \code{hazard = log(2) / median}. Cannot be used together with
#'   \code{e.hazard}.
#' @param e.time Interval boundaries for a piecewise exponential event time,
#'   of the form \code{c(0, t1, ..., Inf)}. Required when \code{e.hazard} or
#'   \code{e.median} is a vector of length > 1. For two groups, supply a list
#'   of two boundary vectors if the change points differ between groups, or a
#'   single vector if they are shared.
#' @param d.hazard Hazard rate(s) for the dropout time. Same structure as
#'   \code{e.hazard}. Set both \code{d.hazard} and \code{d.median} to
#'   \code{NULL} to suppress dropout. Cannot be used together with
#'   \code{d.median}.
#' @param d.median Median dropout time(s). Same structure as \code{e.median}.
#'   Cannot be used together with \code{d.hazard}.
#' @param d.time Interval boundaries for a piecewise exponential dropout time.
#'   Same structure as \code{e.time}.
#' @param seed A single integer for reproducibility. Passed to
#'   \code{dqrng::dqset.seed}. Default is \code{NULL} (no seed).
#'
#' @return A \code{data.frame} with \code{nsim * sum(n)} rows and the
#'   following columns:
#' \describe{
#'   \item{\code{sim}}{Simulation ID (integer, 1 to \code{nsim}).}
#'   \item{\code{group}}{Group label (integer: 1 for control, 2 for
#'     treatment). Always 1 for one-group simulations.}
#'   \item{\code{accrual_time}}{Patient accrual time from study start.}
#'   \item{\code{surv_time}}{Survival time from patient entry.}
#'   \item{\code{dropout_time}}{Dropout time from patient entry
#'     (\code{Inf} when dropout is not modeled).}
#'   \item{\code{tte}}{Observed time-to-event:
#'     \code{pmin(surv_time, dropout_time)}.}
#'   \item{\code{event}}{Event indicator: 1 if \code{surv_time <=
#'     dropout_time}, 0 otherwise.}
#'   \item{\code{calendar_time}}{Calendar time from study start:
#'     \code{accrual_time + tte}.}
#' }
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
#' # Two-group simulation, simple exponential, with dropout
#' df2 <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(100, 100),
#'   a.time   = c(0, 6, 12),
#'   a.rate   = c(8, 12),
#'   e.median = list(18, 24),
#'   d.hazard = list(0.01, 0.01),
#'   seed     = 2
#' )
#' head(df2)
#'
#' # Two-group simulation, piecewise exponential (delayed treatment effect)
#' df3 <- simdata_fast(
#'   nsim     = 100,
#'   n        = c(100, 100),
#'   a.time   = c(0, 12),
#'   a.rate   = 200 / 12,
#'   e.hazard = list(c(0.08, 0.08), c(0.08, 0.04)),
#'   e.time   = c(0, 6, Inf),
#'   seed     = 3
#' )
#' head(df3)
#'
#' # Two-group via total n + allocation ratio
#' df4 <- simdata_fast(
#'   nsim     = 100,
#'   n        = 200,
#'   alloc    = c(1, 1),
#'   a.time   = c(0, 12),
#'   a.rate   = 200 / 12,
#'   e.hazard = list(0.08, 0.05),
#'   seed     = 4
#' )
#' head(df4)
#'
#' @importFrom dqrng dqset.seed dqrexp
#' @export
simdata_fast <- function(nsim     = 1000,
                         n,
                         alloc    = c(1, 1),
                         a.time,
                         a.rate,
                         e.hazard = NULL,
                         e.median = NULL,
                         e.time   = NULL,
                         d.hazard = NULL,
                         d.median = NULL,
                         d.time   = NULL,
                         seed     = NULL) {

  # ------------------------------------------------------------------ #
  #  Seed
  # ------------------------------------------------------------------ #
  if (!is.null(seed)) dqrng::dqset.seed(seed)

  # ------------------------------------------------------------------ #
  #  Resolve group sizes
  # ------------------------------------------------------------------ #
  if (length(n) == 1L) {
    n_grp    <- round(n * alloc / sum(alloc))
    n_groups <- 2L
  } else if (length(n) == 2L) {
    n_grp    <- n
    n_groups <- 2L
  } else {
    stop("'n' must be a scalar (total N) or a vector of length 2 (per-group)")
  }

  is_two_group <- n_groups == 2L && (
    is.list(e.hazard) || is.list(e.median) ||
      is.list(d.hazard) || is.list(d.median)
  )
  if (!is_two_group && length(n) == 2L) {
    stop("Two-group sample sizes supplied but hazard/median parameters are not lists. ",
         "Wrap group-specific parameters in list().")
  }

  if (!is_two_group) {
    n_grp    <- if (length(n) == 1L) n else n[1L]
    n_groups <- 1L
  }

  # ------------------------------------------------------------------ #
  #  Validate and resolve e.hazard / e.median
  # ------------------------------------------------------------------ #
  if (!is.null(e.hazard) && !is.null(e.median)) {
    stop("Specify exactly one of 'e.hazard' and 'e.median', not both")
  }
  if (is.null(e.hazard) && is.null(e.median)) {
    stop("One of 'e.hazard' or 'e.median' must be supplied")
  }
  if (!is.null(e.median)) e.hazard <- convert_median_to_hazard(e.median)

  # ------------------------------------------------------------------ #
  #  Validate and resolve d.hazard / d.median
  # ------------------------------------------------------------------ #
  has_dropout <- !is.null(d.hazard) || !is.null(d.median)
  if (!is.null(d.hazard) && !is.null(d.median)) {
    stop("Specify exactly one of 'd.hazard' and 'd.median', not both")
  }
  if (!is.null(d.median)) d.hazard <- convert_median_to_hazard(d.median)

  # ------------------------------------------------------------------ #
  #  Validate accrual parameters
  # ------------------------------------------------------------------ #
  if (length(a.rate) != length(a.time) - 1L) {
    stop("'a.rate' must have length equal to length(a.time) - 1")
  }
  if (any(a.rate <= 0)) stop("All 'a.rate' values must be positive")
  if (any(diff(a.time) <= 0)) stop("'a.time' must be strictly increasing")

  # Pre-compute accrual cumulative probability table (shared across groups)
  lengths_a <- diff(a.time)
  weights_a <- a.rate * lengths_a
  cum_p_a   <- cumsum(weights_a / sum(weights_a))

  # ------------------------------------------------------------------ #
  #  Simulate each group and combine
  # ------------------------------------------------------------------ #
  if (n_groups == 1L) {
    df <- simulate_group(
      nsim     = nsim,
      n        = n_grp,
      group_id = 1L,
      a.time   = a.time,
      cum_p_a  = cum_p_a,
      e.hazard = if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard,
      e.time   = e.time,
      d.hazard = if (has_dropout) (if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard) else NULL,
      d.time   = d.time
    )
  } else {
    e.time_ctrl <- resolve_time_arg(e.time, 1L)
    e.time_trt  <- resolve_time_arg(e.time, 2L)
    d.time_ctrl <- if (has_dropout) resolve_time_arg(d.time, 1L) else NULL
    d.time_trt  <- if (has_dropout) resolve_time_arg(d.time, 2L) else NULL

    df_ctrl <- simulate_group(
      nsim     = nsim,
      n        = n_grp[1L],
      group_id = 1L,
      a.time   = a.time,
      cum_p_a  = cum_p_a,
      e.hazard = if (is.list(e.hazard)) e.hazard[[1L]] else e.hazard,
      e.time   = e.time_ctrl,
      d.hazard = if (has_dropout) (if (is.list(d.hazard)) d.hazard[[1L]] else d.hazard) else NULL,
      d.time   = d.time_ctrl
    )

    df_trt <- simulate_group(
      nsim     = nsim,
      n        = n_grp[2L],
      group_id = 2L,
      a.time   = a.time,
      cum_p_a  = cum_p_a,
      e.hazard = if (is.list(e.hazard)) e.hazard[[2L]] else e.hazard,
      e.time   = e.time_trt,
      d.hazard = if (has_dropout) (if (is.list(d.hazard)) d.hazard[[2L]] else d.hazard) else NULL,
      d.time   = d.time_trt
    )

    # C++ interleave: replaces rbind() + order()
    df <- interleave_groups(df_ctrl, df_trt, nsim,
                            as.integer(n_grp[1L]), as.integer(n_grp[2L]))
  }

  df
}

# ------------------------------------------------------------------ #
#  Internal helper: simulate one group
# ------------------------------------------------------------------ #
simulate_group <- function(nsim, n, group_id, a.time, cum_p_a,
                           e.hazard, e.time, d.hazard, d.time) {
  total_n <- nsim * n

  # Accrual times via C++ piecewise uniform sampler
  accrual_time <- rpiece_unif_cpp(total_n, a.time, cum_p_a)

  # Survival times
  surv_time <- rpiece_exp_r(total_n, e.hazard, e.time)

  # Dropout times
  dropout_time <- if (!is.null(d.hazard)) {
    rpiece_exp_r(total_n, d.hazard, d.time)
  } else {
    rep(Inf, total_n)
  }

  tte           <- pmin(surv_time, dropout_time)
  event         <- as.integer(surv_time <= dropout_time)
  calendar_time <- accrual_time + tte

  data.frame(
    sim           = rep(seq_len(nsim), each = n),
    group         = group_id,
    accrual_time  = accrual_time,
    surv_time     = surv_time,
    dropout_time  = dropout_time,
    tte           = tte,
    event         = event,
    calendar_time = calendar_time
  )
}

# ------------------------------------------------------------------ #
#  Internal helper: piecewise exponential sampler
# ------------------------------------------------------------------ #
rpiece_exp_r <- function(n, hazard, e.time) {
  if (length(hazard) == 1L) {
    return(dqrng::dqrexp(n, rate = hazard))
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

  # C++ piecewise exponential sampler
  rpiece_exp_cpp(n, hazard, fin_time, cum_haz)
}

# ------------------------------------------------------------------ #
#  Internal helper: resolve e.time / d.time for a specific group
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
    lapply(median_arg, function(m) log(2) / m)
  } else {
    log(2) / median_arg
  }
}
