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

// Forward declarations of the pointer-based implementations (defined below).
void pihe_core_impl(const double*, const int*, const int*, int,
                    std::vector<PiheEvSummary>&, double*);
void pihe_core_strat_impl(const double*, const int*, const int*, const int*,
                          int, std::vector<PiheEvSummary>&, double*);

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

//' Core stratified PiHE hazard ratio computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the quantities of the stratified
//' Pike-Halley Estimator. The at-risk sets are restarted in every stratum and
//' the observed and expected event totals, the score, the information, and
//' the curvature term are summed over strata, which gives the Pike anchor and
//' its Halley correction for the stratified Breslow partial likelihood. Not
//' intended to be called directly by users; use \code{coxph_fast()} with
//' \code{strata} instead.
//'
//' @param time_sorted A numeric vector of pooled follow-up times, sorted by
//'   stratum and in ascending order of time within stratum.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param strata_sorted An integer vector of stratum codes aligned with
//'   \code{time_sorted}; rows of the same stratum must be contiguous.
//'
//' @return A numeric vector of length 4: \code{c(theta_0, U_0, I_0, J_0)}.
//'   Returns a length-4 vector of \code{NA_real_} when the estimate cannot
//'   be computed.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector pihe_core_strat(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    const IntegerVector& strata_sorted
) {
  const int n = time_sorted.size();
  std::vector<PiheEvSummary> ev;
  double out[4];
  pihe_core_strat_impl(time_sorted.begin(), event_sorted.begin(),
                       j_sorted.begin(), strata_sorted.begin(), n, ev, out);
  return NumericVector::create(out[0], out[1], out[2], out[3]);
}

// Shared scan for the unstratified and stratified estimators. When strata is
// a null pointer the whole input is a single stratum. Otherwise rows of the
// same stratum must be contiguous and sorted by time within the stratum; the
// at-risk counts are restarted at every stratum, the observed and expected
// totals and the per-event-time summaries are pooled over strata, and the
// Pike anchor, score, information, and curvature are therefore those of the
// stratified Breslow partial likelihood (a sum over strata).
static void pihe_scan(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    const int* strata,
    int n,
    std::vector<PiheEvSummary>& ev,
    double* out
) {
  out[0] = NA_REAL; out[1] = NA_REAL; out[2] = NA_REAL; out[3] = NA_REAL;
  if (n == 0) return;

  ev.clear();
  if ((int) ev.capacity() < n) ev.reserve(n);

  double O_T = 0.0, O_C = 0.0, E_T = 0.0, E_C = 0.0;

  // Pass 1: one left-to-right scan per stratum block [b, e)
  int b = 0;
  while (b < n) {
    int e = n;
    if (strata != nullptr) {
      const int s = strata[b];
      e = b;
      while (e < n && strata[e] == s) ++e;
    }

    // Initialize at-risk counts for this block
    int n1 = 0, n0 = 0;
    for (int k = b; k < e; ++k) {
      if (j_sorted[k] == 1) ++n1; else ++n0;
    }

    int i = b;
    while (i < e) {
      const double t = time_sorted[i];

      // Consume the tied block at time t
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

    b = e;
  }

  if (O_T == 0.0 || O_C == 0.0 || E_T == 0.0 || E_C == 0.0) return;

  // Pike anchor
  const double theta_0 = (O_T * E_C) / (O_C * E_T);

  // Pass 2: walk saved summaries to compute U_0, I_0, J_0 at theta_0
  double U_0 = 0.0, I_0 = 0.0, J_0 = 0.0;
  const std::size_t K = ev.size();
  for (std::size_t k = 0; k < K; ++k) {
    const PiheEvSummary& es = ev[k];
    const double denom = es.n0k + es.n1k * theta_0;
    if (denom == 0.0) return;

    const double p_k = es.n1k * theta_0 / denom;
    const double q_k = 1.0 - p_k;
    const double pq  = p_k * q_k;
    const double Opq = es.Ok * pq;

    U_0 += es.OTk - es.Ok * p_k;
    I_0 += Opq;
    J_0 += Opq * (1.0 - 2.0 * p_k);
  }

  if (!std::isfinite(I_0) || I_0 == 0.0) return;

  out[0] = theta_0; out[1] = U_0; out[2] = I_0; out[3] = J_0;
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
  pihe_scan(time_sorted, event_sorted, j_sorted, nullptr, n, ev, out);
}

// Stratified counterpart with external linkage. The input must be sorted by
// stratum and by time within stratum (rows of a stratum contiguous).
void pihe_core_strat_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    const int* strata_sorted,
    int n,
    std::vector<PiheEvSummary>& ev,
    double* out
) {
  pihe_scan(time_sorted, event_sorted, j_sorted, strata_sorted, n, ev, out);
}
