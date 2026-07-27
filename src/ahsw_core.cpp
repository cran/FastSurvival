// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
// dH, Rrun, Gfrac are reusable scratch buffers supplied by the caller.
void ahsw_core_impl(const double*, const int*, int, double,
                    std::vector<double>&, std::vector<double>&,
                    std::vector<double>&, double*);

//' Core average hazard with survival weight computation for one group (C++ backend)
//'
//' @description
//' Internal C++ function that computes the average hazard with survival weight
//' (AHSW) of Uno and Horiguchi for a single group, together with the building
//' blocks needed for its asymptotic variance. The average hazard on the window
//' from 0 to \code{tau} is the ratio of the cumulative event probability at
//' \code{tau} to the restricted mean survival time at \code{tau}, both based on
//' the Kaplan-Meier estimate. The function walks a pooled sorted dataset once,
//' accumulating the Nelson-Aalen increment, the Kaplan-Meier survival, the
//' running restricted mean survival time and the at-risk fraction at each
//' distinct event time at or before \code{tau}, and forms the two variance
//' sums used for the ratio (log scale) and the difference (identity scale).
//' Not intended to be called directly by users; use \code{ahsw_fast()} instead.
//'
//' @details
//' Writing F for the cumulative event probability at \code{tau}, R for the
//' restricted mean survival time at \code{tau}, R(t) for the running restricted
//' mean survival time, dH(t) for the Nelson-Aalen increment and G(t) for the
//' at-risk fraction, the average hazard is F / R and the two variance sums are
//' the sum over event times of dH(t) (kernel(t))^2 / G(t). For the average
//' hazard on the log scale the kernel is 1 / F - R(t) / R, matching the
//' \code{v_Q} term of the reference implementation, and for the average hazard
//' on the identity scale the kernel is 1 / R - F R(t) / R^2, matching the
//' \code{v_U} term. The survival function is the right-continuous Kaplan-Meier
//' product, and the running restricted mean survival time integrates the
//' left-hand survival over each grid interval, exactly as in the reference.
//'
//' @param time_sorted A numeric vector of follow-up times sorted in ascending
//'   order.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param tau A numeric truncation time point.
//'
//' @return A numeric vector of length 6: \code{c(F_tau, RMST_tau, AH, v_Q, v_U,
//'   surv_tau)}, where F_tau is the cumulative event probability at \code{tau},
//'   RMST_tau is the restricted mean survival time at \code{tau}, AH is the
//'   average hazard F_tau / RMST_tau, v_Q is the variance sum on the log scale,
//'   v_U is the variance sum on the identity scale, and surv_tau is the
//'   Kaplan-Meier survival at \code{tau}. The variance sums must be divided by
//'   the group size by the caller.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector ahsw_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    double tau
) {
  const int n = time_sorted.size();
  std::vector<double> dH, Rrun, Gfrac;
  double out[6];
  ahsw_core_impl(time_sorted.begin(), event_sorted.begin(), n, tau,
                 dH, Rrun, Gfrac, out);
  return NumericVector::create(out[0], out[1], out[2], out[3], out[4], out[5]);
}

// Pointer-based implementation with external linkage. Writes
// F_tau, RMST_tau, AH, v_Q, v_U, surv_tau into out[0..5]. The three per-event
// scratch vectors are cleared on entry and reused across cells. Algorithm
// identical to the exported wrapper above.
void ahsw_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    int n,
    double tau,
    std::vector<double>& dH,      // Nelson-Aalen increment at each event <= tau
    std::vector<double>& Rrun,    // running RMST to the left of each event time
    std::vector<double>& Gfrac,   // at-risk fraction at each event time
    double* out
) {
  dH.clear();
  Rrun.clear();
  Gfrac.clear();

  double surv      = 1.0;      // right-continuous KM survival
  double rmst      = 0.0;      // running RMST = integral_0^t S
  double prev_time = 0.0;      // previous grid time (starts at 0)
  double surv_prev = 1.0;      // survival on the left of the current interval
  double surv_tau  = 1.0;      // KM survival at tau

  int i = 0;
  while (i < n) {
    const double t = time_sorted[i];
    if (t > tau) break;

    // Consume the tied block at time t
    int d = 0, c = 0;
    int jj = i;
    while (jj < n && time_sorted[jj] == t) {
      ++c;
      if (event_sorted[jj] == 1) ++d;
      ++jj;
    }

    const int nrisk = n - i;   // subjects with time >= t (sorted ascending)

    if (d > 0 && nrisk > 0) {
      // Integrate survival over [prev_time, t] using the left survival
      rmst += (t - prev_time) * surv_prev;

      const double dh = (double)d / (double)nrisk;
      dH.push_back(dh);
      Rrun.push_back(rmst);
      Gfrac.push_back((double)nrisk / (double)n);

      // Update KM survival after this event time
      surv *= (1.0 - dh);
      surv_prev = surv;
      prev_time = t;
    }
    i = jj;
  }

  // Integrate survival over the final stretch up to tau
  if (prev_time < tau) {
    rmst += (tau - prev_time) * surv_prev;
  }
  surv_tau = surv;

  const double F_tau    = 1.0 - surv_tau;
  const double RMST_tau = rmst;
  const double AH       = (RMST_tau > 0.0) ? (F_tau / RMST_tau) : NA_REAL;

  // Second pass: variance sums
  double v_Q = 0.0, v_U = 0.0;
  if (F_tau > 0.0 && RMST_tau > 0.0) {
    const int m = (int)dH.size();
    for (int a = 0; a < m; ++a) {
      if (dH[a] == 0.0 || Gfrac[a] == 0.0) continue;
      const double kerQ = 1.0 / F_tau - Rrun[a] / RMST_tau;
      const double kerU = 1.0 / RMST_tau - F_tau * Rrun[a] / (RMST_tau * RMST_tau);
      v_Q += dH[a] * kerQ * kerQ / Gfrac[a];
      v_U += dH[a] * kerU * kerU / Gfrac[a];
    }
  } else {
    v_Q = NA_REAL;
    v_U = NA_REAL;
  }

  out[0] = F_tau; out[1] = RMST_tau; out[2] = AH;
  out[3] = v_Q;   out[4] = v_U;      out[5] = surv_tau;
}
