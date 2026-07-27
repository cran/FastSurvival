#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Result container for the two-group average hazard ratio core.
struct AHRResult {
  double theta1;
  double theta2;
  double ahr;
  double var_theta1;
  double var_theta2;
  double GL;
  int valid;
};

// Evaluate one group's Kaplan-Meier survival curve and the cumulative Greenwood
// variance of log(S) (scaled by the group size) on a common grid of times.
// The group's observations must be sorted in ascending time order. All group
// observations enter the risk set, so times beyond the grid still contribute to
// the at-risk counts at earlier event times (no left truncation is assumed).
//   gtime, gstatus : right-censored times and 0/1 event indicators for the group
//   grid           : common evaluation times (ascending, unique), right-continuous
//   ng             : group sample size
//   surv, phi      : outputs, length grid.size()
static void ahr_eval_km(const std::vector<double>& gtime,
                        const std::vector<int>& gstatus,
                        const std::vector<double>& grid,
                        int ng,
                        std::vector<double>& surv,
                        std::vector<double>& phi) {
  int n = (int) gtime.size();
  int m = (int) grid.size();
  surv.assign(m, 1.0);
  phi.assign(m, 0.0);

  double s_cur = 1.0;   // running survival probability
  double v_cur = 0.0;   // running cumulative Greenwood variance of log(S)
  int at_risk = n;      // subjects with time >= current event time
  int idx = 0;          // pointer into the (sorted) observations
  int g = 0;            // pointer into the grid

  while (idx < n) {
    double t_cur = gtime[idx];

    // Aggregate ties at t_cur
    int d = 0;      // events at t_cur
    int tied = 0;   // all observations at t_cur (events + censored)
    int j = idx;
    while (j < n && gtime[j] == t_cur) {
      if (gstatus[j] == 1) d++;
      tied++;
      j++;
    }

    // Grid points strictly before t_cur take the pre-update (right-continuous) value
    while (g < m && grid[g] < t_cur) {
      surv[g] = s_cur;
      phi[g] = (double) ng * v_cur;
      g++;
    }

    // Update Kaplan-Meier and Greenwood at the event time
    if (d > 0 && at_risk > 0) {
      s_cur *= 1.0 - (double) d / (double) at_risk;
      if (at_risk - d > 0) {
        v_cur += (double) d / ((double) at_risk * (double) (at_risk - d));
      }
      // if at_risk - d == 0 the Greenwood increment is undefined; it is skipped
    }

    // Grid points at exactly t_cur take the post-update value
    while (g < m && grid[g] == t_cur) {
      surv[g] = s_cur;
      phi[g] = (double) ng * v_cur;
      g++;
    }

    at_risk -= tied;
    idx = j;
  }

  // Grid points after the last observation keep the final value
  while (g < m) {
    surv[g] = s_cur;
    phi[g] = (double) ng * v_cur;
    g++;
  }
}

// Core computation of the Kalbfleisch-Prentice average hazard ratio over [0, tau]
// for two groups, together with the direct (Greenwood-based) variance of theta1.
// Separated from the exported wrapper so it can be reused from other translation
// units (fused use inside simulation loops) without R-level overhead.
static AHRResult ahr_core_impl(const std::vector<double>& time1,
                               const std::vector<int>& status1,
                               const std::vector<double>& time2,
                               const std::vector<int>& status2,
                               const std::vector<double>& grid) {
  int n1 = (int) time1.size();
  int n2 = (int) time2.size();
  int n = n1 + n2;
  double p1 = (double) n1 / (double) n;
  double p2 = (double) n2 / (double) n;
  int m = (int) grid.size();

  std::vector<double> S1, phi1, S2, phi2;
  ahr_eval_km(time1, status1, grid, n1, S1, phi1);
  ahr_eval_km(time2, status2, grid, n2, S2, phi2);

  // Jumps dS_g(t_i) = S_g(t_i) - S_g(t_{i - 1}), with S_g(t_{-1}) = 1
  std::vector<double> dS1(m), dS2(m);
  for (int i = 0; i < m; ++i) {
    double prev1 = (i == 0) ? 1.0 : S1[i - 1];
    double prev2 = (i == 0) ? 1.0 : S2[i - 1];
    dS1[i] = S1[i] - prev1;
    dS2[i] = S2[i] - prev2;
  }

  double GL = S1[m - 1] * S2[m - 1];
  double denom = 1.0 - GL;

  AHRResult res;
  res.GL = GL;

  // No events in [0, tau] (or degenerate): estimate undefined
  if (denom <= 0.0) {
    res.theta1 = res.theta2 = res.ahr = res.var_theta1 = res.var_theta2 = NA_REAL;
    res.valid = 0;
    return res;
  }

  // Point estimate: theta1 = - integral of S2 dS1 over [0, tau], divided by (1 - GL)
  double x1 = 0.0;
  for (int i = 0; i < m; ++i) x1 += S2[i] * dS1[i];

  double theta1 = -x1 / denom;
  double theta2 = 1.0 - theta1;
  res.theta1 = theta1;
  res.theta2 = theta2;

  if (theta1 <= 0.0 || theta2 <= 0.0) {
    res.ahr = NA_REAL;
    res.var_theta1 = NA_REAL;
    res.var_theta2 = NA_REAL;
    res.valid = 0;
    return res;
  }
  res.ahr = theta2 / theta1;

  // Variance pieces (Kalbfleisch-Prentice direct estimator specialised to k = 2)
  double A12 = 0.0, A21 = 0.0;
  for (int i = 0; i < m; ++i) {
    A12 += S1[i] * phi1[i] * dS2[i];
    A21 += S2[i] * phi2[i] * dS1[i];
  }
  A12 *= GL / p1;
  A21 *= GL / p2;

  double C1 = GL * GL * phi1[m - 1] / p1;
  double C2 = GL * GL * phi2[m - 1] / p2;

  // Totals for the tail parts of the double-integral terms
  double T1 = 0.0, T2 = 0.0;
  for (int i = 0; i < m; ++i) {
    T1 += S1[i] * dS2[i];
    T2 += S2[i] * dS1[i];
  }

  // B[1,2,2] and B[2,1,1] via O(m) prefix sums, exploiting the
  // independent-increments structure phi(min(s, u)) of the log-survival
  // covariance. inner(u) = sum_{s <= u} . phi(s) . + phi(u) sum_{s > u} .
  double B122 = 0.0;
  {
    double cum_pref = 0.0, cum_sds = 0.0;
    for (int u = 0; u < m; ++u) {
      cum_pref += S1[u] * phi1[u] * dS2[u];
      cum_sds += S1[u] * dS2[u];
      double inner = cum_pref + phi1[u] * (T1 - cum_sds);
      B122 += S1[u] * inner * dS2[u];
    }
    B122 /= p1;
  }

  double B211 = 0.0;
  {
    double cum_pref = 0.0, cum_sds = 0.0;
    for (int u = 0; u < m; ++u) {
      cum_pref += S2[u] * phi2[u] * dS1[u];
      cum_sds += S2[u] * dS1[u];
      double inner = cum_pref + phi2[u] * (T2 - cum_sds);
      B211 += S2[u] * inner * dS1[u];
    }
    B211 /= p2;
  }

  double Vx11 = B122 + B211 + C1 - 2.0 * A12;
  double VxG1 = A12 - A21 - C1;
  double VG = C1 + C2;

  double Sigma11 = Vx11 + 2.0 * theta1 * VxG1 + theta1 * theta1 * VG;
  Sigma11 /= denom * denom;
  res.var_theta1 = Sigma11 / (double) n;

  // Variance of the second (comparison) group share. The direct estimator does
  // not force Sigma[1,1] == Sigma[2,2] in finite samples, and the AHR package
  // tests each share with its own diagonal entry, so this is computed too.
  double Vx22 = B211 + B122 + C2 - 2.0 * A21;
  double VxG2 = A21 - A12 - C2;
  double Sigma22 = Vx22 + 2.0 * theta2 * VxG2 + theta2 * theta2 * VG;
  Sigma22 /= denom * denom;
  res.var_theta2 = Sigma22 / (double) n;

  res.valid = 1;
  return res;
}

// Pointer-based two-group entry with external linkage. Builds the evaluation
// grid (0, the pooled event times on or before tau, and tau) internally from
// the two groups' sorted observations and writes
//   theta1, theta2, ahr, var_theta1, var_theta2, GL, valid
// into out[0..6]. The two groups' observations must each be sorted in ascending
// time order. Provided so a fused simulation loop can call the average hazard
// ratio core directly on reusable group-split buffers, matching the grid that
// the exported wrapper builds in R.
void ahr_core_two_impl(const double* t0, const int* d0, int n0,
                       const double* t1, const int* d1, int n1,
                       double tau, double* out) {
  std::vector<double> grid;
  grid.reserve((size_t) (n0 + n1 + 2));
  grid.push_back(0.0);
  for (int i = 0; i < n0; ++i) if (d0[i] == 1 && t0[i] <= tau) grid.push_back(t0[i]);
  for (int i = 0; i < n1; ++i) if (d1[i] == 1 && t1[i] <= tau) grid.push_back(t1[i]);
  grid.push_back(tau);
  std::sort(grid.begin(), grid.end());
  grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

  std::vector<double> tt0(t0, t0 + n0), tt1(t1, t1 + n1);
  std::vector<int>    dd0(d0, d0 + n0), dd1(d1, d1 + n1);

  AHRResult r = ahr_core_impl(tt0, dd0, tt1, dd1, grid);
  out[0] = r.theta1;
  out[1] = r.theta2;
  out[2] = r.ahr;
  out[3] = r.var_theta1;
  out[4] = r.var_theta2;
  out[5] = r.GL;
  out[6] = (double) r.valid;
}

//' Two-group average hazard ratio core (internal)
//'
//' Internal Rcpp core for \code{\link{ahr_fast}}. Computes the
//' Kalbfleisch-Prentice average hazard ratio over the interval from 0 to
//' \code{tau} for two groups, with the direct Greenwood-based variance of the
//' first group share \code{theta1}. The two groups' observations must each be
//' sorted in ascending time order and \code{grid} must be sorted and unique.
//'
//' @param time1 sorted right-censored times for the reference group
//' @param status1 0/1 event indicators for the reference group
//' @param time2 sorted right-censored times for the comparison group
//' @param status2 0/1 event indicators for the comparison group
//' @param grid common evaluation times (sorted, unique, including 0 and tau)
//' @return A named numeric vector with elements theta1, theta2, ahr,
//'   var.theta1, var.theta2, GL and valid.
//' @keywords internal
// [[Rcpp::export]]
NumericVector ahr_core(NumericVector time1, IntegerVector status1,
                       NumericVector time2, IntegerVector status2,
                       NumericVector grid) {
  std::vector<double> t1(time1.begin(), time1.end());
  std::vector<int> d1(status1.begin(), status1.end());
  std::vector<double> t2(time2.begin(), time2.end());
  std::vector<int> d2(status2.begin(), status2.end());
  std::vector<double> gr(grid.begin(), grid.end());

  AHRResult r = ahr_core_impl(t1, d1, t2, d2, gr);

  return NumericVector::create(
    _["theta1"] = r.theta1,
    _["theta2"] = r.theta2,
    _["ahr"] = r.ahr,
    _["var.theta1"] = r.var_theta1,
    _["var.theta2"] = r.var_theta2,
    _["GL"] = r.GL,
    _["valid"] = (double) r.valid
  );
}
