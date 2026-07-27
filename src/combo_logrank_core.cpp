// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Forward declaration of the pointer-based implementation (defined below).
// u_out has length nw, v_out has length nw*nw (row-major), w is a reusable
// scratch buffer of length nw. Returns O1.
double combo_logrank_core_impl(const double*, const int*, const int*, int,
                               const double*, const double*, int,
                               double*, double*, std::vector<double>&);

//' Core max-combo weighted log-rank computation on pooled sorted vectors (C++ backend)
//'
//' @description
//' Internal C++ function that computes a vector of Fleming-Harrington weighted
//' log-rank numerators and their full covariance matrix in a single pass over a
//' pooled sorted dataset. Each weight is a Fleming-Harrington G(rho, gamma)
//' weight evaluated from the left-continuous pooled Kaplan-Meier estimate, with
//' the pairs supplied through \code{rho_vec} and \code{gamma_vec}. The variance
//' increment at each event time does not depend on the weight, so it is computed
//' once and reused for every numerator and every pair of weights. The diagonal
//' of the returned covariance matrix holds the per-scheme variances, and the
//' off-diagonal entries hold the between-scheme covariances used to form the
//' correlation matrix for the max-combo p-value. Tied event times are processed
//' atomically. Not intended to be called directly by users; use the max-combo
//' R wrapper instead.
//'
//' @details
//' For a weight indexed by k the numerator is \code{U_k = sum w_k (O1 - E1)} and
//' the variance is \code{V_kk = sum w_k^2 v}, where \code{v} is the
//' hypergeometric variance increment \code{d n1 (n - n1) (n - d) /
//' (n^2 (n - 1))}. The between-scheme covariance is \code{V_ab = sum w_a w_b v}.
//' Here \code{S_minus} is the left-continuous pooled Kaplan-Meier estimate just
//' prior to each event time, initialized at 1, and the weight is
//' \code{S_minus^rho (1 - S_minus)^gamma}. The convention is the same as the
//' single-scheme weighted core: the numerator is negative when the treatment
//' group is favored.
//'
//' @param time_sorted A numeric vector of pooled follow-up times sorted in
//'   ascending order.
//' @param event_sorted An integer vector of event indicators (1 = event,
//'   0 = censored), aligned with \code{time_sorted}.
//' @param j_sorted An integer vector of group indicators (1 = treatment,
//'   0 = control), aligned with \code{time_sorted}.
//' @param rho_vec A numeric vector of Fleming-Harrington first parameters, one
//'   per weight scheme.
//' @param gamma_vec A numeric vector of Fleming-Harrington second parameters,
//'   one per weight scheme, aligned with \code{rho_vec}.
//'
//' @return A list with three elements: \code{O1}, the observed number of events
//'   in the treatment group; \code{U}, a numeric vector of weighted log-rank
//'   numerators (one per weight scheme); and \code{V}, the numeric covariance
//'   matrix of the numerators.
//'
//' @keywords internal
// [[Rcpp::export]]
List combo_logrank_core(
    const NumericVector& time_sorted,
    const IntegerVector& event_sorted,
    const IntegerVector& j_sorted,
    const NumericVector& rho_vec,
    const NumericVector& gamma_vec
) {
  const int n = time_sorted.size();
  const int nw = rho_vec.size();

  NumericVector u_acc(nw);
  NumericMatrix v_acc(nw, nw);
  std::vector<double> w(nw);

  // v_acc is column-major (R matrix); the impl writes row-major into a flat
  // buffer, but since V is symmetric the layouts coincide, so we can pass
  // v_acc.begin() directly.
  const double O1 = combo_logrank_core_impl(
    time_sorted.begin(), event_sorted.begin(), j_sorted.begin(), n,
    rho_vec.begin(), gamma_vec.begin(), nw,
    u_acc.begin(), v_acc.begin(), w);

  return List::create(
    _["O1"] = O1,
    _["U"]  = u_acc,
    _["V"]  = v_acc
  );
}

// Pointer-based implementation with external linkage. Writes the nw weighted
// numerators into u_out[0..nw-1] and the nw*nw covariance matrix into
// v_out (symmetric, so row-major and column-major coincide). w is a reusable
// scratch buffer of length nw. Returns O1. The caller must zero u_out and
// v_out before the call. Algorithm identical to the exported wrapper above.
double combo_logrank_core_impl(
    const double* time_sorted,
    const int* event_sorted,
    const int* j_sorted,
    int n,
    const double* rho_vec,
    const double* gamma_vec,
    int nw,
    double* u_out,
    double* v_out,
    std::vector<double>& w
) {
  // At-risk counts per group, initialized to the group totals
  int n1_init = 0, n0_init = 0;
  for (int k = 0; k < n; ++k) {
    if (j_sorted[k] == 1) ++n1_init; else ++n0_init;
  }

  double O1 = 0.0;

  int n1 = n1_init, n0 = n0_init;
  double s_minus = 1.0;   // left-continuous pooled KM just prior to t
  double s       = 1.0;   // right-continuous pooled KM

  int i = 0;
  while (i < n) {
    const double t = time_sorted[i];

    // Consume tied block at time t, splitting deaths and censorings per group
    int d1 = 0, d0 = 0, c1 = 0, c0 = 0;
    int jj = i;
    while (jj < n && time_sorted[jj] == t) {
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

      const double e1 = dd * dn1 / dnj;
      const double v1 = dd * dn1 * (dnj - dn1) * (dnj - dd) /
                        (dnj * dnj * (dnj - 1.0));
      const double oe = (double)d1 - e1;

      // Fleming-Harrington weight per scheme from the left-continuous KM
      for (int a = 0; a < nw; ++a) {
        w[a] = std::pow(s_minus, rho_vec[a]) *
               std::pow(1.0 - s_minus, gamma_vec[a]);
        u_out[a] += w[a] * oe;
      }

      // Symmetric covariance increments w_a w_b v1
      for (int a = 0; a < nw; ++a) {
        const double wa_v = w[a] * v1;
        for (int b = a; b < nw; ++b) {
          const double inc = wa_v * w[b];
          v_out[a * nw + b] += inc;
          if (b != a) v_out[b * nw + a] += inc;
        }
      }

      O1 += d1;
    }

    // Update KM after using s_minus at this time: s_minus -> s -> next s_minus
    if (d > 0 && nj > 0) {
      s *= (1.0 - (double)d / (double)nj);
    }
    s_minus = s;

    // Decrement at-risk counts by this block's size
    n1 -= c1;
    n0 -= c0;
    i   = jj;
  }

  return O1;
}
