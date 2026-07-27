#' Fast Log-Rank Test for Two-Group Survival Data
#'
#' @description
#' Computes the log-rank test statistic for comparing survival curves between
#' two groups. Returns either a one-sided Z-score or a two-sided chi-square
#' statistic. The C++ backend uses a two-pointer merge scan over pooled sorted
#' vectors, eliminating the \code{rank()} call that dominates the pure-R
#' implementation. When a \code{strata} argument is supplied, the stratified
#' log-rank test is computed instead, matching
#' \code{\link[survival]{survdiff}} with a \code{strata()} term. A non-default
#' \code{weight} argument selects a weighted log-rank test (Fleming-Harrington,
#' modestly-weighted, Gehan-Breslow, or Tarone-Ware) for non-proportional
#' hazards.
#'
#' @details
#' The log-rank statistic is computed as:
#'
#' Z = (O_1 - E_1) / sqrt(V_1)
#'
#' where O_1 is the observed number of events in the treatment group, E_1 is
#' the expected number under the null hypothesis of equal survival, and V_1 is
#' the hypergeometric variance. Tied event times are handled correctly: all
#' subjects sharing the same event time form a tied block, and the block is
#' processed atomically in the two-pointer merge.
#'
#' When \code{strata} is \code{NULL} (default), the ordinary two-group log-rank
#' test is computed by the C++ core \code{logrank_core}, which walks the pooled
#' sorted data with a single two-pointer scan, maintaining running at-risk
#' counts per group. No rank vector is constructed, so the dominant O(n log n)
#' cost of \code{rank()} in the pure-R version is removed.
#'
#' When \code{strata} is supplied, the stratified log-rank test is computed by
#' the C++ core \code{stratified_logrank_core}. The contributions O_1, E_1, and
#' V_1 are accumulated within each stratum and then summed across strata, so
#' the overall statistic is Z = (sum O_1 - sum E_1) / sqrt(sum V_1). A stratum
#' that contains only one group contributes zero to all three totals, the same
#' convention used by \code{\link[survival]{survdiff}}. This is the standard
#' stratified log-rank test, equivalent to a log-rank test that conditions on
#' the stratum at each event time.
#'
#' When a non-default \code{weight} is requested, the weighted log-rank test is
#' computed by the C++ core \code{weighted_logrank_core}. The statistic is
#' Z = U / sqrt(V) with U = sum w_j (O_1j - E_1j) and
#' V = sum w_j^2 n_0j n_1j O_j (n_j - O_j) / (n_j^2 (n_j - 1)), where the weight
#' w_j is one of the schemes named by \code{weight}. The Fleming-Harrington and
#' modestly-weighted schemes use the left-continuous pooled Kaplan-Meier
#' estimate S(t-), initialized at 1 and updated after each event time, so the
#' first event always has S(t-) = 1. The modestly-weighted scheme determines
#' its weight cap in a first pass over the event times before accumulating U
#' and V in a second pass; the other schemes accumulate in a single pass. When
#' both \code{weight} and \code{strata} are supplied, the stratified weighted
#' log-rank test is computed by the C++ core
#' \code{stratified_weighted_logrank_core}: each stratum is an independent
#' weighted log-rank test whose weights come from that stratum's own pooled
#' Kaplan-Meier estimate, and the per-stratum U and V are summed before
#' standardizing once as Z = sum U / sqrt(sum V).
#'
#' When \code{presorted = TRUE}, the input vectors are assumed to be sorted and
#' the internal \code{order()} call is skipped. In the unstratified case the
#' assumed order is ascending \code{time}. In the stratified case the assumed
#' order is by stratum first and by ascending \code{time} within each stratum,
#' so that rows of the same stratum are contiguous. When \code{presorted =
#' FALSE} (default), sorting is handled internally. In simulation loops where
#' the data are generated in the required order, setting \code{presorted =
#' TRUE} avoids one O(n log n) pass.
#'
#' The returned object has class \code{"survdiff_fast"} and is a length-one
#' numeric value (Z-score or chi-square) with the underlying counts O_0, E_0,
#' O_1, E_1, V_1, the requested \code{side}, and the total sample size stored
#' as attributes. When \code{strata} is supplied, the number of strata is also
#' stored as the attribute \code{strata}. A \code{print()} method formats the
#' result similarly to \code{print(survival::survdiff(...))}, displaying
#' observed and expected event counts for both the control and treatment
#' groups.
#'
#' @param time A numeric vector of follow-up times for all subjects.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param side An integer, either 1 or 2 (default 2). If \code{side = 1}, returns
#'   the standardized log-rank statistic (Z-score), defined as
#'   \code{(O_1 - E_1) / sqrt(V_1)} for the treatment group, so the Z-score is
#'   negative when the treatment group has fewer events than expected (a
#'   protective treatment effect). If \code{side = 2}, returns the chi-square
#'   statistic (Z^2).
#' @param presorted A logical value. If \code{TRUE}, \code{time}, \code{event},
#'   and \code{group} (and \code{strata} when supplied) are assumed to be
#'   already sorted in the required order, and the internal \code{order()} call
#'   is skipped. If \code{FALSE} (default), sorting is handled internally. See
#'   Details for the required order in the stratified case.
#' @param strata An optional vector of stratum labels aligned with \code{time}.
#'   If \code{NULL} (default), the ordinary log-rank test is computed and the
#'   behavior is identical to earlier versions of this function. If supplied,
#'   the stratified log-rank test is computed, matching
#'   \code{\link[survival]{survdiff}} with a \code{strata()} term. Any type that
#'   supports equality comparison is accepted. May be combined with a non-default
#'   \code{weight} to obtain a stratified weighted log-rank test.
#' @param weight A character string naming the weight scheme. \code{"logrank"}
#'   (default) is the ordinary unweighted log-rank test and reproduces the
#'   behavior of earlier versions of this function exactly. \code{"fh"} is the
#'   Fleming-Harrington G(rho, gamma) test with weight
#'   \code{S(t-)^rho (1 - S(t-))^gamma}. \code{"mwlrt"} is the modestly-weighted
#'   log-rank test of Magirr and Burman with weight
#'   \code{1 / max(S(t-), S(t_star))}. \code{"gehan"} is the Gehan-Breslow test
#'   with weight equal to the at-risk count, and \code{"tarone-ware"} uses the
#'   square root of the at-risk count. Here \code{S(t-)} is the left-continuous
#'   pooled Kaplan-Meier estimate just prior to each event time.
#' @param rho A numeric Fleming-Harrington first parameter, used only when
#'   \code{weight = "fh"}. Defaults to 0.
#' @param gamma A numeric Fleming-Harrington second parameter, used only when
#'   \code{weight = "fh"}. Defaults to 0. The pair \code{rho = 0, gamma = 0}
#'   reproduces the ordinary log-rank test, and \code{rho = 0, gamma = 1} is the
#'   Fleming-Harrington G(0, 1) test for delayed effects.
#' @param t_star A single non-negative numeric value, the timepoint of the
#'   modestly-weighted log-rank test. Required only when \code{weight =
#'   "mwlrt"}. The weight is capped at \code{1 / S(t_star)}, where
#'   \code{S(t_star)} is the pooled Kaplan-Meier value at \code{t_star}.
#'   A value of 0 yields the ordinary log-rank test.
#'
#' @return An object of class \code{"survdiff_fast"}, which is a length-one
#'   numeric value with attributes \code{O0}, \code{E0}, \code{O1}, \code{E1},
#'   \code{V1}, \code{side}, and \code{n}, plus \code{strata} (the number of
#'   strata) when \code{strata} is supplied, or \code{weight} (the scheme name)
#'   when a non-default \code{weight} is used. The numeric value is the Z-score
#'   when \code{side = 1}, or the chi-square statistic when \code{side = 2}.
#'   For a weighted test the value is U / sqrt(V) (or its square), \code{V1}
#'   holds the weighted variance V, \code{O0} and \code{O1} hold the raw
#'   observed event counts, and \code{E0} and \code{E1} are \code{NA} because a
#'   single unweighted expected count is not defined for a weighted test.
#'   Returns \code{NA_real_} (still with class \code{"survdiff_fast"}) when the
#'   variance is zero (e.g., all events in one group).
#'
#' @examples
#' library(survival)
#'
#' # Two-sided test: compare with survdiff
#' fit <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 2)
#' fit
#'
#' chisq_ref <- survdiff(Surv(futime, fustat) ~ rx, data = ovarian)$chisq
#' cat("survdiff_fast chi-square:", as.numeric(fit), "\n")
#' cat("survdiff      chi-square:", chisq_ref,       "\n")
#'
#' # One-sided test (Z-score)
#' survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1)
#'
#' # presorted = TRUE: sort once outside, reuse inside a loop
#' ord <- order(ovarian$futime)
#' survdiff_fast(ovarian$futime[ord], ovarian$fustat[ord], ovarian$rx[ord],
#'               control = 2, side = 2, presorted = TRUE)
#'
#' # Stratified log-rank test: compare with survdiff + strata()
#' fit_str <- survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2,
#'                          side = 2, strata = ovarian$resid.ds)
#' chisq_str <- survdiff(Surv(futime, fustat) ~ rx + strata(resid.ds),
#'                       data = ovarian)$chisq
#' cat("stratified survdiff_fast:", as.numeric(fit_str), "\n")
#' cat("stratified survdiff     :", chisq_str,           "\n")
#'
#' # Weighted log-rank tests for non-proportional hazards
#' # Fleming-Harrington G(0, 1), emphasizing late differences
#' survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1,
#'               weight = "fh", rho = 0, gamma = 1)
#'
#' # Modestly-weighted log-rank test with t_star = 365
#' survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1,
#'               weight = "mwlrt", t_star = 365)
#'
#' # Stratified weighted log-rank test: Fleming-Harrington G(0,1) within strata
#' survdiff_fast(ovarian$futime, ovarian$fustat, ovarian$rx, 2, side = 1,
#'               weight = "fh", rho = 0, gamma = 1, strata = ovarian$resid.ds)
#'
#' \donttest{
#' # Speed comparison against survdiff()
#' if (requireNamespace("microbenchmark", quietly = TRUE)) {
#'   microbenchmark::microbenchmark(
#'     survdiff_fast = survdiff_fast(ovarian$futime, ovarian$fustat,
#'                                   ovarian$rx, 2, side = 2),
#'     survdiff      = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
#'     times = 1000
#'   )
#' }
#' }
#'
#' @seealso
#' \code{\link[survival]{survdiff}} for the standard implementation.
#' \code{\link{print.survdiff_fast}} for the print method.
#'
#' @references
#' Gehan, E. A. (1965). A generalized Wilcoxon test for comparing arbitrarily
#' single-censored samples. \emph{Biometrika}, \emph{52}, 203-223.
#'
#' Mantel, N. (1966). Evaluation of survival data and two new rank order
#' statistics arising in its consideration. \emph{Cancer Chemotherapy
#' Reports}, \emph{50}(3), 163-170.
#'
#' Tarone, R. E., & Ware, J. (1977). On distribution-free tests for equality of
#' survival distributions. \emph{Biometrika}, \emph{64}, 156-160.
#'
#' Fleming, T. R., & Harrington, D. P. (1991). \emph{Counting Processes and
#' Survival Analysis}. New York: John Wiley & Sons.
#'
#' Magirr, D., & Burman, C.-F. (2019). Modestly weighted logrank tests.
#' \emph{Statistics in Medicine}, \emph{38}(20), 3782-3790.
#'
#' @export
survdiff_fast <- function(time, event, group, control, side = 2,
                          presorted = FALSE, strata = NULL,
                          weight = c("logrank", "fh", "mwlrt", "gehan",
                                     "tarone-ware"),
                          rho = 0, gamma = 0, t_star = NULL) {
  weight <- match.arg(weight)

  # Input validation
  if (length(time) != length(event) || length(time) != length(group)) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  if (!side %in% c(1L, 2L)) {
    stop("'side' must be either 1 (one-sided) or 2 (two-sided)")
  }
  if (sum(event) == 0L) {
    stop("No events observed in the data")
  }
  if (weight == "mwlrt" &&
      (is.null(t_star) || length(t_star) != 1L || !is.finite(t_star) ||
       t_star < 0)) {
    stop("'t_star' must be a single non-negative value when weight = 'mwlrt'")
  }

  n_total <- length(time)

  # ------------------------------------------------------------------ #
  #  Weighted log-rank path (new). Taken only when a non-default weight
  #  is requested; the unweighted paths below are unchanged.
  # ------------------------------------------------------------------ #
  if (weight != "logrank") {
    scheme <- switch(weight,
                     fh          = 0L,
                     mwlrt       = 1L,
                     gehan       = 2L,
                     `tarone-ware` = 3L)
    t_star_v <- if (weight == "mwlrt") t_star else 0

    # Treatment indicator: 1 = treatment, 0 = control
    j <- as.integer(group != control)

    use_strata <- !is.null(strata)
    if (use_strata) {
      if (length(strata) != n_total) {
        stop("'strata' must have the same length as 'time'")
      }
      # Map stratum labels to contiguous integers 1..S
      strata_int <- match(strata, sort(unique(strata)))
      # Sort by stratum first, then by time within stratum
      if (!presorted) {
        ord        <- order(strata_int, time)
        time       <- time[ord]
        event      <- as.integer(event[ord])
        j          <- j[ord]
        strata_int <- strata_int[ord]
      } else {
        event <- as.integer(event)
      }
      # C++ core: per-stratum weighted scan, totals summed -> c(O1, U, V)
      res <- stratified_weighted_logrank_core(time, event, j, strata_int,
                                              scheme, rho, gamma, t_star_v)
      n_str <- length(unique(strata_int))
    } else {
      if (!presorted) {
        ord   <- order(time)
        time  <- time[ord]
        event <- as.integer(event[ord])
        j     <- j[ord]
      } else {
        event <- as.integer(event)
      }
      # C++ core: single scan (two passes for mwlrt) -> c(O1, U, V)
      res <- weighted_logrank_core(time, event, j, scheme, rho, gamma, t_star_v)
      n_str <- NULL
    }

    O1 <- res[1L]
    U  <- res[2L]
    V  <- res[3L]

    total_events <- sum(event)
    O0           <- total_events - O1

    base_attr <- list(
      O0 = O0, E0 = NA_real_, O1 = O1, E1 = NA_real_, V1 = V,
      side = side, n = n_total, weight = weight, control = control
    )
    if (use_strata) base_attr$strata <- n_str

    if (!is.finite(V) || V == 0) {
      return(do.call(structure, c(list(NA_real_), base_attr,
                                  list(class = "survdiff_fast"))))
    }

    WLR <- U / sqrt(V)
    val <- if (side == 1L) WLR else WLR^2

    return(do.call(structure, c(list(val), base_attr,
                                list(class = "survdiff_fast"))))
  }

  # ------------------------------------------------------------------ #
  #  Stratified log-rank path (new). Taken only when 'strata' is given;
  #  the unstratified path below is unchanged and identical to earlier
  #  versions of this function.
  # ------------------------------------------------------------------ #
  if (!is.null(strata)) {
    if (length(strata) != n_total) {
      stop("'strata' must have the same length as 'time'")
    }

    # Treatment indicator: 1 = treatment, 0 = control
    j <- as.integer(group != control)

    # Map stratum labels to contiguous integers 1..S
    strata_int <- match(strata, sort(unique(strata)))

    # Sort by stratum first, then by time within stratum, when not presorted
    if (!presorted) {
      ord        <- order(strata_int, time)
      time       <- time[ord]
      event      <- as.integer(event[ord])
      j          <- j[ord]
      strata_int <- strata_int[ord]
    } else {
      event <- as.integer(event)
    }

    # C++ core: per-stratum scan, totals summed across strata -> c(O1, E1, V1)
    ovr <- stratified_logrank_core(time, event, j, strata_int)
    O1  <- ovr[1L]
    E1  <- ovr[2L]
    V1  <- ovr[3L]

    total_events <- sum(event)
    O0           <- total_events - O1
    E0           <- total_events - E1
    n_str        <- length(unique(strata_int))

    if (!is.finite(V1) || V1 == 0) {
      return(structure(
        NA_real_,
        O0 = O0, E0 = E0, O1 = O1, E1 = E1, V1 = V1,
        side = side, n = n_total, strata = n_str, control = control,
        class = "survdiff_fast"
      ))
    }

    LR  <- (O1 - E1) / sqrt(V1)
    val <- if (side == 1L) LR else LR^2

    return(structure(
      val,
      O0 = O0, E0 = E0, O1 = O1, E1 = E1, V1 = V1,
      side = side, n = n_total, strata = n_str, control = control,
      class = "survdiff_fast"
    ))
  }

  # ------------------------------------------------------------------ #
  #  Unstratified log-rank path (unchanged)
  # ------------------------------------------------------------------ #

  # Treatment indicator: 1 = treatment, 0 = control
  j <- as.integer(group != control)

  # Sort pooled data by time when not presorted
  if (!presorted) {
    ord   <- order(time)
    time  <- time[ord]
    event <- as.integer(event[ord])
    j     <- j[ord]
  } else {
    event <- as.integer(event)
  }

  # C++ core: single pass over pooled sorted data -> c(O1, E1, V1)
  ovr <- logrank_core(time, event, j)
  O1  <- ovr[1L]
  E1  <- ovr[2L]
  V1  <- ovr[3L]

  # Control-group counts: O0 = total events - O1, E0 = total events - E1
  total_events <- sum(event)
  O0           <- total_events - O1
  E0           <- total_events - E1

  if (!is.finite(V1) || V1 == 0) {
    return(structure(
      NA_real_,
      O0 = O0, E0 = E0, O1 = O1, E1 = E1, V1 = V1,
      side = side, n = n_total, control = control,
      class = "survdiff_fast"
    ))
  }

  LR  <- (O1 - E1) / sqrt(V1)
  val <- if (side == 1L) LR else LR^2

  structure(
    val,
    O0 = O0, E0 = E0, O1 = O1, E1 = E1, V1 = V1,
    side = side, n = n_total, control = control,
    class = "survdiff_fast"
  )
}
