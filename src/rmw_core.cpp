// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
void rmw_core_impl(const double*, const int*, const int*, int, double,
                   double*);

//' Core robust modestly-weighted log-rank computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes, in a single pass over a pooled sorted
//' dataset, the two component statistics of the robust modestly-weighted (rMW)
//' log-rank test of Magirr and Ohrn together with their null covariance. The
//' first component is the standard log-rank statistic (weight one at every
//' event time); the second is a modestly-weighted log-rank statistic with
//' weight \code{min(1 / S(t-), 1 / s_star)}, where \code{S(t-)} is the
//' left-continuous pooled Kaplan-Meier estimate just prior to each event time
//' and \code{s_star} is a survival-probability threshold. The rMW test
//' rejects for an extreme value of the maximum of the two standardized
//' components, so the covariance of the two numerators under the null is
//' required to recover the joint distribution. Tied event times are processed
//' atomically. Not intended to be called directly by users; use
//' \code{rmw_fast()} instead.
//'
//' @details
//' The standard log-rank numerator is \code{U_lr = sum (d1 - e1)} with variance
//' \code{V_lr = sum var_d}, where \code{d1} is the number of treatment-group
//' events at an event time, \code{e1 = d n1 / n} is the expected count, and
//' \code{var_d = d n1 n0 (n - d) / (n^2 (n - 1))} is the hypergeometric
//' variance of \code{d1}. The modestly-weighted numerator is
//' \code{U_mw = sum w (d1 - e1)} with variance \code{V_mw = sum w^2 var_d}.
//' Because the standard log-rank weight is one, the null covariance of the two
//' numerators reduces to \code{C = sum w var_d}. The modestly-weighted weight
//' is capped at \code{1 / s_star}, a constant, so unlike the timepoint
//' parameterization used by \code{weighted_logrank_core} no first pass is
//' needed and all quantities accumulate in one scan. When \code{s_star} is not
//' positive the cap is treated as infinite, giving the uncapped
//' \code{1 / S(t-)} weight. Setting \code{s_star = 1} caps the weight at one,
//' so the modestly-weighted component equals the standard log-rank component
//' and the two are perfectly correlated.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param s_star A numeric survival-probability threshold for the
//'   modestly-weighted component. The weight is capped at \code{1 / s_star},
//'   so a value of \code{0.5} caps the weight at 2.
//'
//' @return A numeric vector of length 6:
//'   \code{c(O1, U_lr, V_lr, U_mw, V_mw, C)}, where O1 is the observed number
//'   of events in the treatment group, U_lr and V_lr are the standard log-rank
//'   numerator and variance, U_mw and V_mw are the modestly-weighted numerator
//'   and variance, and C is the null covariance of the two numerators.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector rmw_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    double s_star
) {
  const int n = time_sorted.size();
  double out[6];
  rmw_core_impl(time_sorted.begin(), event_sorted.begin(), j_sorted.begin(),
                n, s_star, out);
  return NumericVector::create(out[0], out[1], out[2], out[3], out[4], out[5]);
}

// Pointer-based implementation with external linkage. Writes
// O1, U_lr, V_lr, U_mw, V_mw, C into out[0..5]. The algorithm is identical to
// the exported wrapper above and is provided so a fused simulation loop can
// call it directly on reusable buffers.
void rmw_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    int n,
    double s_star,
    double* out
) {
  // Weight cap 1 / s_star. A non-positive s_star means no cap.
  const double cap = (s_star > 0.0) ? (1.0 / s_star) : R_PosInf;

  // At-risk counts per group, initialized to the group totals
  int n1 = 0, n0 = 0;
  for (int k = 0; k < n; ++k) {
    if (j_sorted[k] == 1) ++n1; else ++n0;
  }

  double O1 = 0.0;
  double U_lr = 0.0, V_lr = 0.0;
  double U_mw = 0.0, V_mw = 0.0;
  double C = 0.0;

  double s_minus = 1.0;   // left-continuous pooled KM just prior to t
  double s_km    = 1.0;   // right-continuous pooled KM

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

      const double e1    = dd * dn1 / dnj;
      const double var_d = dd * dn1 * dn0 * (dnj - dd) /
                           (dnj * dnj * (dnj - 1.0));
      const double dev   = (double)d1 - e1;

      // Modestly-weighted weight from the left-continuous KM, capped at 1/s_star
      double w = 1.0 / s_minus;
      if (w > cap) w = cap;

      O1   += d1;
      U_lr += dev;
      V_lr += var_d;
      U_mw += w * dev;
      V_mw += w * w * var_d;
      C    += w * var_d;
    }

    // Update KM after using s_minus at this time: s_minus -> s_km -> next s_minus
    if (d > 0 && nj > 0) {
      s_km *= (1.0 - (double)d / (double)nj);
    }
    s_minus = s_km;

    // Decrement at-risk counts by this block's size
    n1 -= c1;
    n0 -= c0;
    i   = j;
  }

  out[0] = O1;
  out[1] = U_lr; out[2] = V_lr;
  out[3] = U_mw; out[4] = V_mw;
  out[5] = C;
}
