#' Fast Log-Rank Test for Two-Group Survival Data
#'
#' @description
#' Computes the log-rank test statistic for comparing survival curves between
#' two groups. Returns either a one-sided Z-score or a two-sided chi-square
#' statistic. The C++ backend uses a two-pointer merge scan over pooled sorted
#' vectors, eliminating the \code{rank()} call that dominates the pure-R
#' implementation.
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
#' When \code{presorted = TRUE}, the input vectors are assumed to be sorted in
#' ascending order of \code{time} and the internal \code{order()} call is
#' skipped. When \code{presorted = FALSE} (default), sorting is handled
#' internally. In simulation loops where the data are generated in sorted order,
#' setting \code{presorted = TRUE} avoids one O(n log n) pass.
#'
#' The C++ core (\code{logrank_core}) walks the pooled sorted data with a
#' single two-pointer scan, maintaining running at-risk counts per group. No
#' rank vector is constructed, so the dominant O(n log n) cost of \code{rank()}
#' in the pure-R version is removed.
#'
#' The returned object has class \code{"survdiff_fast"} and is a length-one
#' numeric value (Z-score or chi-square) with the underlying counts O_0, E_0,
#' O_1, E_1, V_1, the requested \code{side}, and the total sample size stored
#' as attributes. A \code{print()} method formats the result similarly to
#' \code{print(survival::survdiff(...))}, displaying observed and expected
#' event counts for both the control and treatment groups.
#'
#' @param time A numeric vector of follow-up times for all subjects.
#' @param event An integer or numeric vector of event indicators
#'   (1 = event, 0 = censored), aligned with \code{time}.
#' @param group A vector of group labels aligned with \code{time}.
#' @param control A scalar value indicating which level of \code{group}
#'   represents the control group.
#' @param side An integer, either 1 or 2. If \code{side = 1}, returns the
#'   standardized log-rank statistic (Z-score, positive when the treatment
#'   group has fewer events than expected). If \code{side = 2}, returns the
#'   chi-square statistic (Z^2).
#' @param presorted A logical value. If \code{TRUE}, \code{time}, \code{event},
#'   and \code{group} are assumed to be sorted in ascending order of
#'   \code{time}, and the internal \code{order()} call is skipped. If
#'   \code{FALSE} (default), sorting is handled internally.
#'
#' @return An object of class \code{"survdiff_fast"}, which is a length-one
#'   numeric value with attributes \code{O0}, \code{E0}, \code{O1}, \code{E1},
#'   \code{V1}, \code{side}, and \code{n}. The numeric value is the Z-score
#'   when \code{side = 1}, or the chi-square statistic when \code{side = 2}.
#'   Returns \code{NA_real_} (still with class \code{"survdiff_fast"}) when
#'   \code{V1 = 0} (e.g., all events in one group).
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
#' \donttest{
#' library(microbenchmark)
#' microbenchmark(
#'   survdiff_fast = survdiff_fast(ovarian$futime, ovarian$fustat,
#'                                 ovarian$rx, 2, side = 2),
#'   survdiff      = survdiff(Surv(futime, fustat) ~ rx, data = ovarian),
#'   times = 1000
#' )
#' }
#'
#' @seealso
#' \code{\link[survival]{survdiff}} for the standard implementation.
#' \code{\link{print.survdiff_fast}} for the print method.
#'
#' @references
#' Mantel, N. (1966). Evaluation of survival data and two new rank order
#' statistics arising in its consideration. \emph{Cancer Chemotherapy
#' Reports}, \emph{50}(3), 163-170.
#'
#' Peto, R., & Peto, J. (1972). Asymptotically efficient rank invariant test
#' procedures. \emph{Journal of the Royal Statistical Society. Series A
#' (General)}, \emph{135}(2), 185-198.
#'
#' @export
survdiff_fast <- function(time, event, group, control, side,
                          presorted = FALSE) {
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

  n_total <- length(time)

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
      side = side, n = n_total,
      class = "survdiff_fast"
    ))
  }

  LR  <- (O1 - E1) / sqrt(V1)
  val <- if (side == 1L) LR else LR^2

  structure(
    val,
    O0 = O0, E0 = E0, O1 = O1, E1 = E1, V1 = V1,
    side = side, n = n_total,
    class = "survdiff_fast"
  )
}
