// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
void stratified_weighted_logrank_core_impl(const double*, const int*,
                                           const int*, const int*, int,
                                           int, double, double, double,
                                           double*);

//' Core stratified weighted log-rank computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the stratified weighted log-rank totals
//' from a dataset sorted by stratum first and by time within stratum. Each
//' stratum is processed as an independent weighted log-rank test using the same
//' logic as \code{weighted_logrank_core}, with weights derived from the
//' within-stratum pooled Kaplan-Meier estimate, and the per-stratum numerator U
//' and variance V are accumulated into the returned totals. The overall
//' statistic is Z = sum U / sqrt(sum V). This matches the stratified weighted
//' log-rank test of Magirr, which sums the per-stratum U and V and standardizes
//' once. Not intended to be called directly by users; use
//' \code{survdiff_fast()} with both \code{weight} and \code{strata} instead.
//'
//' @details
//' The scheme codes match \code{weighted_logrank_core}: 0 = Fleming-Harrington
//' G(rho, gamma); 1 = modestly-weighted with cap \code{1 / S(t_star)} computed
//' within each stratum; 2 = Gehan-Breslow; 3 = Tarone-Ware. The left-continuous
//' pooled Kaplan-Meier estimate is restarted at 1 at the beginning of every
//' stratum, so the weights of each stratum depend only on that stratum's
//' pooled data. For the modestly-weighted scheme the weight cap is determined
//' in a first pass within each stratum before accumulation.
//'
//' @param time_sorted A numeric vector of follow-up times sorted by stratum
//'   first and by ascending time within each stratum.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param strata_sorted An integer vector of stratum identifiers, aligned with
//'   \code{time_sorted}. Rows of the same stratum must be contiguous.
//' @param scheme An integer weight-scheme code: 0 = Fleming-Harrington,
//'   1 = modestly-weighted, 2 = Gehan-Breslow, 3 = Tarone-Ware.
//' @param rho A numeric Fleming-Harrington first parameter (scheme 0).
//' @param gamma A numeric Fleming-Harrington second parameter (scheme 0).
//' @param t_star A numeric timepoint for the modestly-weighted scheme
//'   (scheme 1).
//'
//' @return A numeric vector of length 3: \code{c(O1, U, V)}, where O1 is the
//'   pooled observed number of events in the treatment group, U is the summed
//'   weighted numerator, and V is the summed weighted variance.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector stratified_weighted_logrank_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    const IntegerVector& strata_sorted,
    int scheme,
    double rho,
    double gamma,
    double t_star
) {
  const int n = time_sorted.size();
  double out[3];
  stratified_weighted_logrank_core_impl(
    time_sorted.begin(), event_sorted.begin(), j_sorted.begin(),
    strata_sorted.begin(), n, scheme, rho, gamma, t_star, out);
  return NumericVector::create(out[0], out[1], out[2]);
}

// Pointer-based implementation with external linkage. Writes O1, U, V into
// out[0..2]. Algorithm identical to the exported wrapper above.
void stratified_weighted_logrank_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    const int* strata_sorted,
    int n,
    int scheme,
    double rho,
    double gamma,
    double t_star,
    double* out
) {
  double O1_tot = 0.0, U_tot = 0.0, V_tot = 0.0;

  // Walk the data stratum block by stratum block.
  int b = 0;
  while (b < n) {
    const int s = strata_sorted[b];

    // Find the end of the current stratum block [b, e)
    int e = b;
    while (e < n && strata_sorted[e] == s) ++e;

    // Count at-risk per group in this stratum
    int n1_init = 0, n0_init = 0;
    for (int k = b; k < e; ++k) {
      if (j_sorted[k] == 1) ++n1_init; else ++n0_init;
    }

    // ---- First pass (modestly-weighted only): within-stratum weight cap ---
    double max_weight = 1.0;
    if (scheme == 1 && t_star > 0.0) {
      int nrisk = n1_init + n0_init;
      double s_km = 1.0;
      double s_star = 1.0;   // within-stratum pooled KM at t_star: product over events <= t_star
      int i = b;
      while (i < e) {
        const double t = time_sorted[i];
        int d = 0, c = 0;
        int jj = i;
        while (jj < e && time_sorted[jj] == t) {
          ++c;
          if (event_sorted[jj] == 1) ++d;
          ++jj;
        }
        if (d > 0 && nrisk > 0) {
          s_km *= (1.0 - (double)d / (double)nrisk);
        }
        if (t < t_star) s_star = s_km;   // within-stratum KM just before t_star (events strictly < t_star)
        nrisk -= c;
        i = jj;
      }
      if (s_star > 0.0) {
        max_weight = 1.0 / s_star;
      }
    }

    // ---- Main accumulation pass within the stratum ------------------------
    int n1 = n1_init, n0 = n0_init;
    double s_minus = 1.0;   // left-continuous pooled KM within stratum
    double s_km    = 1.0;   // right-continuous pooled KM within stratum

    int i = b;
    while (i < e) {
      const double t = time_sorted[i];

      int d1 = 0, d0 = 0, c1 = 0, c0 = 0;
      int jj = i;
      while (jj < e && time_sorted[jj] == t) {
        if (j_sorted[jj] == 1) {
          ++c1;
          if (event_sorted[jj] == 1) ++d1;
        } else {
          ++c0;
          if (event_sorted[jj] == 1) ++d0;
        }
        ++jj;
      }

      const int d  = d1 + d0;
      const int nj = n1 + n0;

      if (d > 0 && nj > 1) {
        const double dn1 = (double)n1;
        const double dnj = (double)nj;
        const double dd  = (double)d;

        double w;
        if (scheme == 0) {
          w = std::pow(s_minus, rho) * std::pow(1.0 - s_minus, gamma);
        } else if (scheme == 1) {
          w = 1.0 / s_minus;
          if (w > max_weight) w = max_weight;
        } else if (scheme == 2) {
          w = dnj;
        } else {
          w = std::sqrt(dnj);
        }

        const double e1 = dd * dn1 / dnj;
        const double v1 = dd * dn1 * (dnj - dn1) * (dnj - dd) /
                          (dnj * dnj * (dnj - 1.0));

        O1_tot += d1;
        U_tot  += w * ((double)d1 - e1);
        V_tot  += w * w * v1;
      }

      if (d > 0 && nj > 0) {
        s_km *= (1.0 - (double)d / (double)nj);
      }
      s_minus = s_km;

      n1 -= c1;
      n0 -= c0;
      i   = jj;
    }

    b = e;
  }

  out[0] = O1_tot; out[1] = U_tot; out[2] = V_tot;
}
