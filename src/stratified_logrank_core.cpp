// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
void stratified_logrank_core_impl(const double*, const int*, const int*,
                                  const int*, int, double*);

//' Core stratified log-rank computation on stratum-blocked sorted vectors
//'
//' @description
//' Internal C++ function that computes the stratified log-rank totals O1, E1,
//' and V1 from a dataset sorted by stratum first and by time within stratum.
//' Each stratum is processed as an independent contiguous block using the
//' same single-pass two-group scan as \code{logrank_core}, and the per-stratum
//' O1, E1, and V1 are accumulated into the returned totals. Tied event times
//' are processed atomically within each stratum. A stratum that contains only
//' one group contributes zero to all three totals, matching the behavior of
//' \code{\link[survival]{survdiff}} with a \code{strata()} term. Not intended
//' to be called directly by users; use \code{survdiff_fast()} with the
//' \code{strata} argument instead.
//'
//' @param time_sorted A numeric vector of follow-up times sorted by stratum
//'   first and by ascending time within each stratum.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param strata_sorted An integer vector of stratum identifiers, aligned with
//'   \code{time_sorted}. Rows of the same stratum must be contiguous, which is
//'   guaranteed when the data are sorted by stratum first.
//'
//' @return A numeric vector of length 3: \code{c(O1, E1, V1)}, the stratified
//'   observed events, expected events, and variance for the treatment group.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector stratified_logrank_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    const IntegerVector& strata_sorted
) {
  const int n = time_sorted.size();
  double out[3];
  stratified_logrank_core_impl(time_sorted.begin(), event_sorted.begin(),
                               j_sorted.begin(), strata_sorted.begin(), n, out);
  return NumericVector::create(out[0], out[1], out[2]);
}

// Pointer-based implementation with external linkage. Writes O1, E1, V1 into
// out[0..2]. Algorithm identical to the exported wrapper above.
void stratified_logrank_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    const int* strata_sorted,
    int n,
    double* out
) {
  double O1 = 0.0, E1 = 0.0, V1 = 0.0;

  // Walk the data stratum block by stratum block. Rows of the same stratum
  // are contiguous, so a change in strata_sorted marks a block boundary.
  int b = 0;
  while (b < n) {
    const int s = strata_sorted[b];

    // Find the end of the current stratum block [b, e)
    int e = b;
    while (e < n && strata_sorted[e] == s) ++e;

    // Initialize at-risk counters for this stratum
    int n1 = 0, n0 = 0;
    for (int k = b; k < e; ++k) {
      if (j_sorted[k] == 1) ++n1; else ++n0;
    }

    // Single-pass log-rank scan within the stratum block
    int i = b;
    while (i < e) {
      const double t = time_sorted[i];

      // Consume tied block at time t, splitting deaths and censorings per group
      int d1 = 0, d0 = 0, c1 = 0, c0 = 0;
      int j = i;
      while (j < e && time_sorted[j] == t) {
        if (j_sorted[j] == 1) {
          ++c1;
          if (event_sorted[j] == 1) ++d1;
        } else {
          ++c0;
          if (event_sorted[j] == 1) ++d0;
        }
        ++j;
      }

      const int d  = d1 + d0;
      const int nj = n1 + n0;

      if (d > 0 && nj > 1) {
        const double dn1 = (double)n1;
        const double dn0 = (double)n0;
        const double dnj = (double)nj;
        const double dd  = (double)d;

        O1 += d1;
        E1 += dd * dn1 / dnj;
        V1 += dd * dn1 * dn0 * (dnj - dd) / (dnj * dnj * (dnj - 1.0));
      }

      // Decrement at-risk counts by this block's size
      n1 -= c1;
      n0 -= c0;
      i   = j;
    }

    b = e;
  }

  out[0] = O1; out[1] = E1; out[2] = V1;
}
