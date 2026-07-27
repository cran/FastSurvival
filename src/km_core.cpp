// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
void km_core_impl(const double*, const double*, int, double, double*);

//' Core Kaplan-Meier computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the Kaplan-Meier survival estimate and
//' the Greenwood variance sum at a specified time point. Binary search locates
//' the cutoff index; a single left-to-right scan groups tied times by distinct
//' time, applying one Kaplan-Meier factor and one Greenwood increment per
//' distinct event time. Grouping the ties follows the risk-set convention of
//' \code{survival::survfit} and makes the result independent of the order of
//' events and censorings within a tied time. Accepts the event vector as
//' a numeric (double) vector to avoid an \code{as.integer()} copy in R.
//' Uses \code{std::log1p(-d/n_risk)} for the KM log-product, which is more
//' accurate than \code{std::log(1 - d/n_risk)} near a full drop and avoids one
//' subtraction per distinct event time. Not intended to be called directly by
//' users; use \code{survfit_fast()} instead.
//'
//' @param t_sorted A numeric vector of follow-up times sorted in ascending
//'   order.
//' @param e_sorted A numeric vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{t_sorted}. Accepted as numeric to
//'   avoid an integer-conversion copy on the R side.
//' @param t_eval A single numeric value specifying the evaluation time point.
//'
//' @return A numeric vector of length 2: \code{c(surv, gw_sum)}, where
//'   \code{surv} is the KM estimate and \code{gw_sum} is the Greenwood
//'   variance sum used to compute SE[S(t)] = surv * sqrt(gw_sum).
//'   Returns \code{c(1.0, 0.0)} when no events occur up to \code{t_eval}.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector km_core(
    const NumericVector& t_sorted,
    const NumericVector& e_sorted,
    double t_eval
) {
  const int n = t_sorted.size();
  double out[2];
  km_core_impl(t_sorted.begin(), e_sorted.begin(), n, t_eval, out);
  return NumericVector::create(out[0], out[1]);
}

// Pointer-based implementation with external linkage. Writes surv, gw_sum into
// out[0..1]. Algorithm identical to the exported wrapper above.
void km_core_impl(
    const double* t_sorted,
    const double* e_sorted,
    int n,
    double t_eval,
    double* out
) {
  if (n == 0) { out[0] = R_NaReal; out[1] = 0.0; return; }

  // Binary search: last index m such that t_sorted[m-1] <= t_eval
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (t_sorted[mid] <= t_eval) lo = mid + 1; else hi = mid;
  }
  const int m = lo;

  if (m == 0) { out[0] = 1.0; out[1] = 0.0; return; }

  // Single scan grouping tied times: at each distinct time, count the events d
  // and the total number of records, apply one (1 - d / n_risk) factor, then
  // remove the whole tie group from the risk set. This follows the Kaplan-Meier
  // risk-set convention of survival::survfit and is independent of the order of
  // events and censorings within a tied time. n_risk decreases from n along the
  // scan; it is maintained as a running counter.
  double log_surv = 0.0;
  double gw_sum   = 0.0;
  bool   any_event = false;
  double n_risk = (double)n;

  int i = 0;
  while (i < m) {
    const double t_i = t_sorted[i];
    int d   = 0;
    int cnt = 0;
    while (i < m && t_sorted[i] == t_i) {
      if (e_sorted[i] >= 0.5) ++d;
      ++cnt;
      ++i;
    }
    if (d > 0) {
      any_event = true;
      const double frac = (double)d / n_risk;
      log_surv += std::log1p(-frac);
      gw_sum   += frac / (n_risk - (double)d);
    }
    n_risk -= (double)cnt;
  }

  if (!any_event) { out[0] = 1.0; out[1] = 0.0; return; }

  out[0] = std::exp(log_surv); out[1] = gw_sum;
}
