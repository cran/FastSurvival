#' Build a Set of Two-Group Trial Scenarios for Design Exploration
#'
#' @description
#' Assembles one or more two-group time-to-event scenarios into a single
#' \code{scenario_fast} object for design-stage exploration. Each scenario is a
#' complete set of \code{\link{simdata_fast}} arguments, so the same object can
#' be drawn with \code{\link{plot.scenario_fast}}, summarized with
#' \code{\link{print.scenario_fast}}, and passed to \code{\link{simdata_fast}}
#' to generate data. This is a design-helper that only collects and merges
#' arguments; it performs no simulation and no compiled computation.
#'
#' @details
#' Each element of \code{scenarios} is a named list of arguments that defines one
#' scenario. Those arguments override the shared arguments in \code{shared}, so
#' that parameters held constant across scenarios (sample size, accrual, dropout)
#' are written once in \code{shared} and only the varying parameters are written
#' per scenario. The survival specification follows \code{\link{simdata_fast}}
#' exactly: \code{e.hazard} or \code{e.median} given as a two-element list
#' (control first, treatment second) for the two groups, with \code{e.time}
#' supplying the breakpoints of a piecewise-exponential hazard (last element
#' \code{Inf}). A single hazard or median gives an exponential group.
#'
#' Two optional fields may appear inside a scenario list and are interpreted by
#' this function rather than passed on: \code{label} sets the scenario label, and
#' \code{null} (logical) flags a null scenario for later filtering. Every other
#' field is treated as a \code{\link{simdata_fast}} argument. The scenario label
#' is taken from \code{labels}, then from the \code{label} field, then from the
#' names of \code{scenarios}, then from a default \code{"Scenario k"}.
#'
#' A factorial set of scenarios is built by constructing the \code{scenarios}
#' list with the base tools, for example mapping over the rows of an
#' \code{expand.grid} of the parameters that vary.
#'
#' @param scenarios A non-empty list. Each element is a named list of
#'   \code{\link{simdata_fast}} arguments (optionally with \code{label} and
#'   \code{null}) that defines one scenario.
#' @param shared A named list of \code{\link{simdata_fast}} arguments common to
#'   every scenario. Scenario-specific arguments take precedence.
#' @param labels Optional character vector of scenario labels, one per scenario.
#'
#' @return An object of class \code{scenario_fast}: a list with element
#'   \code{scenarios} (a named list, one entry per scenario, each holding the
#'   merged argument list \code{args}, the \code{label}, and the \code{null}
#'   flag) and element \code{shared} (the shared arguments).
#'
#' @examples
#' # Three two-group scenarios sharing accrual and sample size: a proportional-
#' # hazards case, a delayed-effect case, and a crossing-hazards case.
#' scn <- gen_scenario_fast(
#'   scenarios = list(
#'     "Proportional"   = list(e.median = list(12, 18)),
#'     "Delayed effect" = list(
#'       e.hazard = list(log(2) / 12, c(log(2) / 12, log(2) / 22)),
#'       e.time   = c(0, 6, Inf)
#'     ),
#'     "Crossing"       = list(
#'       e.hazard = list(log(2) / 12, c(log(2) / 7, log(2) / 24)),
#'       e.time   = c(0, 5, Inf)
#'     )
#'   ),
#'   shared = list(n = c(150, 150), a.time = c(0, 12), a.rate = 300 / 12)
#' )
#'
#' print(scn)
#' plot(scn)
#'
#' # Generate data for every scenario with a short composition step. From here
#' # analysis_fast() and simsummary_fast() are applied per scenario as usual.
#' sim_list <- lapply(scn$scenarios, function(s) {
#'   do.call(simdata_fast, c(s$args, list(nsim = 100, seed = 1)))
#' })
#' nrow(sim_list[[1]])
#'
#' @seealso \code{\link{simdata_fast}}, \code{\link{plot.scenario_fast}},
#'   \code{\link{print.scenario_fast}}
#'
#' @export
gen_scenario_fast <- function(scenarios, shared = list(), labels = NULL) {
  if (missing(scenarios) || !is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list of per-scenario argument lists")
  }
  if (!is.list(shared)) {
    stop("'shared' must be a named list of shared arguments")
  }
  n_scn     <- length(scenarios)
  scn_names <- names(scenarios)
  if (is.null(scn_names)) scn_names <- rep("", n_scn)

  if (!is.null(labels)) {
    if (length(labels) != n_scn) {
      stop("length(labels) must equal length(scenarios)")
    }
    base_lab <- as.character(labels)
  } else {
    base_lab <- ifelse(nzchar(scn_names), scn_names,
                       paste("Scenario", seq_len(n_scn)))
  }

  out <- vector("list", n_scn)
  for (k in seq_len(n_scn)) {
    sk <- scenarios[[k]]
    if (!is.list(sk)) {
      stop(sprintf("Scenario %d must be a named list of arguments", k))
    }
    null_flag    <- isTRUE(sk[["null"]])
    sk[["null"]] <- NULL
    lab_k <- if (!is.null(sk[["label"]])) as.character(sk[["label"]]) else base_lab[k]
    sk[["label"]] <- NULL

    merged <- shared
    for (nm in names(sk)) merged[[nm]] <- sk[[nm]]

    if (is.null(merged[["e.hazard"]]) && is.null(merged[["e.median"]])) {
      stop(sprintf("Scenario '%s' has no 'e.hazard' or 'e.median'", lab_k))
    }
    out[[k]] <- list(args = merged, label = lab_k, null = null_flag)
  }
  names(out) <- vapply(out, function(s) s$label, character(1L))

  structure(list(scenarios = out, shared = shared), class = "scenario_fast")
}

# ------------------------------------------------------------------ #
#  Internal helpers: analytic piecewise-exponential survival and HR
# ------------------------------------------------------------------ #

# Piecewise-constant hazard evaluated at 'times'. 'haz' is a length-K hazard
# vector and 'etime' is the length-(K + 1) breakpoint vector with etime[1] == 0
# and the last element Inf, matching simdata_fast. A single hazard is constant.
pwe_hazard <- function(times, haz, etime) {
  if (length(haz) == 1L) return(rep(as.numeric(haz), length(times)))
  k_int <- length(haz)
  brk   <- etime[seq_len(k_int)]
  idx   <- findInterval(times, brk)
  idx[idx < 1L] <- 1L
  haz[idx]
}

# Cumulative hazard at 'times' for the same parametrization. Built so that the
# break-point cumulative values equal simdata_fast's piecewise_precompute().
pwe_cumhaz <- function(times, haz, etime) {
  if (length(haz) == 1L) return(as.numeric(haz) * pmax(times, 0))
  k_int       <- length(haz)
  brk         <- etime[seq_len(k_int)]
  seglen      <- diff(brk)
  cumh_at_brk <- c(0, cumsum(haz[seq_len(k_int - 1L)] * seglen))
  idx         <- findInterval(times, brk)
  idx[idx < 1L] <- 1L
  cumh_at_brk[idx] + haz[idx] * (times - brk[idx])
}

# Resolve one group's hazard vector and breakpoints from a scenario's args.
# group_idx is 1 for control and 2 for treatment.
resolve_group_haz <- function(spec_haz, spec_med, etime, group_idx) {
  if (!is.null(spec_haz)) {
    haz <- if (is.list(spec_haz)) spec_haz[[group_idx]] else spec_haz
  } else {
    med <- if (is.list(spec_med)) spec_med[[group_idx]] else spec_med
    haz <- log(2) / med
  }
  if (is.list(haz)) {
    stop("scenario plotting does not support per-subgroup survival specs")
  }
  haz <- as.numeric(haz)
  et  <- if (is.list(etime)) etime[[group_idx]] else etime
  if (length(haz) > 1L && is.null(et)) {
    stop("'e.time' is required for a piecewise (vector) hazard or median spec")
  }
  list(haz = haz, etime = et)
}

# Evaluate the analytic control and treatment survival curves and the hazard
# ratio over a time grid for one scenario. 'two' is FALSE for a single-arm spec.
scenario_eval <- function(args, times) {
  spec_haz <- args[["e.hazard"]]
  spec_med <- args[["e.median"]]
  etime    <- args[["e.time"]]

  two <- (is.list(spec_haz) && length(spec_haz) >= 2L) ||
         (is.list(spec_med) && length(spec_med) >= 2L)

  g_c <- resolve_group_haz(spec_haz, spec_med, etime, 1L)
  s_c <- exp(-pwe_cumhaz(times, g_c$haz, g_c$etime))

  if (two) {
    g_t <- resolve_group_haz(spec_haz, spec_med, etime, 2L)
    s_t <- exp(-pwe_cumhaz(times, g_t$haz, g_t$etime))
    h_c <- pwe_hazard(times, g_c$haz, g_c$etime)
    h_t <- pwe_hazard(times, g_t$haz, g_t$etime)
    hr  <- h_t / h_c
  } else {
    s_t <- NULL
    hr  <- NULL
  }
  list(times = times, s_c = s_c, s_t = s_t, hr = hr, two = two)
}

# A reasonable default time horizon, derived from breakpoints, medians, hazards,
# and accrual windows across all scenarios.
default_tmax <- function(scn) {
  cand <- numeric(0)
  for (s in scn) {
    a   <- s$args
    et  <- a[["e.time"]]
    ets <- if (is.list(et)) unlist(et) else et
    ets <- ets[is.finite(ets)]
    if (length(ets)) cand <- c(cand, max(ets) * 1.5)
    med  <- a[["e.median"]]
    meds <- if (is.list(med)) unlist(med) else med
    if (length(meds)) cand <- c(cand, max(meds) * 2)
    haz  <- a[["e.hazard"]]
    hazs <- if (is.list(haz)) unlist(haz) else haz
    hazs <- hazs[hazs > 0]
    if (length(hazs)) cand <- c(cand, (log(2) / min(hazs)) * 2)
    at  <- a[["a.time"]]
    atf <- at[is.finite(at)]
    if (length(atf)) cand <- c(cand, max(atf))
  }
  if (length(cand) == 0L) return(24)
  max(cand)
}

# Median read off a survival curve by linear interpolation. Returns NA when the
# curve does not reach 0.5 within the grid.
curve_median <- function(times, surv) {
  if (is.null(surv) || all(surv > 0.5)) return(NA_real_)
  i <- which(surv <= 0.5)[1L]
  if (i == 1L) return(times[1L])
  t0 <- times[i - 1L]; t1 <- times[i]
  s0 <- surv[i - 1L];  s1 <- surv[i]
  t0 + (0.5 - s0) * (t1 - t0) / (s1 - s0)
}
