// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
void logrank_core_impl(const double*, const int*, const int*, int, double*);

//' Core log-rank computation on pooled sorted vectors (C++ backend)
//'
//' @description
//' Internal C++ function that computes O1, E1, and V1 for the log-rank test
//' directly from a pooled sorted dataset and a group indicator. Eliminates
//' the R-level group-splitting copies (\code{time[is1]}, \code{time[!is1]},
//' \code{event[is1]}, \code{event[!is1]}) of the previous design by walking
//' the pooled vector once and updating per-group at-risk counters in C++.
//' Tied event times are processed atomically. Not intended to be called
//' directly by users; use \code{survdiff_fast()} instead.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//'
//' @return A numeric vector of length 3: \code{c(O1, E1, V1)}.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector logrank_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted
) {
  const int n = time_sorted.size();
  double out[3];
  logrank_core_impl(time_sorted.begin(), event_sorted.begin(),
                    j_sorted.begin(), n, out);
  return NumericVector::create(out[0], out[1], out[2]);
}

// Pointer-based implementation with external linkage. Writes O1, E1, V1 into
// out[0..2]. The fused analysis loop calls this directly on its reusable
// buffers to avoid allocating an Rcpp vector per (sim, look, population) cell.
// The algorithm is identical to the exported wrapper above.
void logrank_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    int n,
    double* out
) {
  // Initialize at-risk counters: count totals per group from j_sorted
  int n1 = 0, n0 = 0;
  for (int k = 0; k < n; ++k) {
    if (j_sorted[k] == 1) ++n1; else ++n0;
  }

  double O1 = 0.0, E1 = 0.0, V1 = 0.0;

  int i = 0;
  while (i < n) {
    const double t = time_sorted[i];

    // Consume tied block at time t, splitting deaths and censorings per group
    int d1 = 0, d0 = 0, c1 = 0, c0 = 0;
    int j = i;
    while (j < n && time_sorted[j] == t) {
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

  out[0] = O1; out[1] = E1; out[2] = V1;
}
