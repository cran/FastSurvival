// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Per-event-time summary for the PiHE two-pass computation. Defined at file
// scope so the pointer-based implementation can take a reusable buffer of it
// from the caller (the fused analysis loop) and avoid reallocating per cell.
struct PiheEvSummary {
  double n1k;
  double n0k;
  double OTk;
  double Ok;
};

// Forward declaration of the pointer-based implementation (defined below).
void pihe_core_impl(const double*, const int*, const int*, int,
                    std::vector<PiheEvSummary>&, double*);

//' Core PiHE hazard ratio computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the Pike-Halley Estimator (PiHE) for
//' the hazard ratio. Accepts pooled sorted vectors plus an integer group
//' indicator, performs group splitting and the two-pointer merge scan
//' entirely in C++, and returns the quantities needed for the Halley
//' correction and Wald interval. Not intended to be called directly by
//' users; use \code{coxph_fast()} instead.
//'
//' Implementation notes: Pass 1 accumulates pooled scalars
//' (O_T, O_C, E_T, E_C) and stores per-distinct-event-time summaries in a
//' single struct array (better cache locality than four parallel vectors)
//' with a single \code{reserve(n)} call (no event-count prepass over the
//' input). Pass 2 walks the saved summaries once at the Pike anchor
//' theta_0 to compute U_0, I_0, and J_0.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//'
//' @return A numeric vector of length 4: \code{c(theta_0, U_0, I_0, J_0)}.
//'   Returns a length-4 vector of \code{NA_real_} when the estimate cannot
//'   be computed.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector pihe_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted
) {
  const int n = time_sorted.size();
  std::vector<PiheEvSummary> ev;
  double out[4];
  pihe_core_impl(time_sorted.begin(), event_sorted.begin(),
                 j_sorted.begin(), n, ev, out);
  return NumericVector::create(out[0], out[1], out[2], out[3]);
}

// Pointer-based implementation with external linkage. Writes
// theta_0, U_0, I_0, J_0 into out[0..3], or four NA values when the estimate
// cannot be computed. The per-event-time summaries are stored in the
// caller-supplied buffer ev, which is cleared on entry and reused across cells
// to avoid a reserve(n) allocation per call. Algorithm identical to the
// exported wrapper above.
void pihe_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    int n,
    std::vector<PiheEvSummary>& ev,
    double* out
) {
  out[0] = NA_REAL; out[1] = NA_REAL; out[2] = NA_REAL; out[3] = NA_REAL;
  if (n == 0) return;

  // Initialize at-risk counts: count totals per group in one short loop
  int n1 = 0, n0 = 0;
  for (int k = 0; k < n; ++k) {
    if (j_sorted[k] == 1) ++n1; else ++n0;
  }

  ev.clear();
  if ((int) ev.capacity() < n) ev.reserve(n);

  double O_T = 0.0, O_C = 0.0, E_T = 0.0, E_C = 0.0;

  // Pass 1: single left-to-right scan over pooled sorted data
  int i = 0;
  while (i < n) {
    const double t = time_sorted[i];

    // Consume the tied block at time t
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

    if (d > 0 && nj > 0) {
      const double dn1 = (double)n1;
      const double dn0 = (double)n0;
      const double dnj = (double)nj;
      const double dd  = (double)d;
      const double dd1 = (double)d1;

      O_T += dd1;
      O_C += dd - dd1;
      E_T += dd * dn1 / dnj;
      E_C += dd * dn0 / dnj;

      ev.push_back(PiheEvSummary{dn1, dn0, dd1, dd});
    }

    // Decrement at-risk counts by block size
    n1 -= c1;
    n0 -= c0;
    i   = j;
  }

  if (O_T == 0.0 || O_C == 0.0 || E_T == 0.0 || E_C == 0.0) return;

  // Pike anchor
  const double theta_0 = (O_T * E_C) / (O_C * E_T);

  // Pass 2: walk saved summaries to compute U_0, I_0, J_0 at theta_0
  double U_0 = 0.0, I_0 = 0.0, J_0 = 0.0;
  const std::size_t K = ev.size();
  for (std::size_t k = 0; k < K; ++k) {
    const PiheEvSummary& e = ev[k];
    const double denom = e.n0k + e.n1k * theta_0;
    if (denom == 0.0) return;

    const double p_k = e.n1k * theta_0 / denom;
    const double q_k = 1.0 - p_k;
    const double pq  = p_k * q_k;
    const double Opq = e.Ok * pq;

    U_0 += e.OTk - e.Ok * p_k;
    I_0 += Opq;
    J_0 += Opq * (1.0 - 2.0 * p_k);
  }

  if (!std::isfinite(I_0) || I_0 == 0.0) return;

  out[0] = theta_0; out[1] = U_0; out[2] = I_0; out[3] = J_0;
}
