// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Core Kaplan-Meier computation (C++ backend)
//'
//' @description
//' Internal C++ function that computes the Kaplan-Meier survival estimate and
//' the Greenwood variance sum at a specified time point. Binary search locates
//' the cutoff index; a single left-to-right scan accumulates the KM product
//' and Greenwood sum over event positions only. Accepts the event vector as
//' a numeric (double) vector to avoid an \code{as.integer()} copy in R.
//' Uses \code{std::log1p(-1/n_risk)} for the KM log-product, which is more
//' accurate than \code{std::log(1 - 1/n_risk)} near n_risk = 1 and avoids one
//' subtraction per event. Not intended to be called directly by users; use
//' \code{survfit_fast()} instead.
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
  if (n == 0) return NumericVector::create(R_NaReal, 0.0);

  // Binary search: last index m such that t_sorted[m-1] <= t_eval
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (t_sorted[mid] <= t_eval) lo = mid + 1; else hi = mid;
  }
  const int m = lo;

  if (m == 0) return NumericVector::create(1.0, 0.0);

  // Single scan: accumulate log(KM) and Greenwood sum at event positions only.
  // n_risk decreases from n down to (n - m + 1) along the scan; we maintain
  // it as a running counter rather than recomputing (double)(n - i) each step.
  double log_surv = 0.0;
  double gw_sum   = 0.0;
  bool   any_event = false;
  double n_risk = (double)n;

  for (int i = 0; i < m; ++i) {
    if (e_sorted[i] >= 0.5) {
      any_event = true;
      const double inv = 1.0 / n_risk;
      log_surv += std::log1p(-inv);
      gw_sum   += inv / (n_risk - 1.0);
    }
    n_risk -= 1.0;
  }

  if (!any_event) return NumericVector::create(1.0, 0.0);

  return NumericVector::create(std::exp(log_surv), gw_sum);
}
