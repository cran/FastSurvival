// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(dqrng)]]
#include <Rcpp.h>
#include <dqrng.h>
#include <cmath>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

// ------------------------------------------------------------------ //
//  Internal: piecewise uniform sampler (accrual times)
// ------------------------------------------------------------------ //
// Equivalent to R's rpiece_unif() but avoids two dqrunif() calls,
// findInterval(), pmax/pmin, and vector index operations by performing
// all steps in a single C++ loop.
//
// Strategy: draw one uniform u per subject. Use u to select the interval
// (via the pre-computed cumulative probability table), then rescale the
// fractional part of u within that interval to get the position.
// This reuses a single draw for both interval selection and within-interval
// positioning, cutting the number of RNG calls by half.
//
// [[Rcpp::export]]
NumericVector rpiece_unif_cpp(
    int n,
    const NumericVector& a_time,
    const NumericVector& cum_p    // cumulative selection probabilities (length = n_int)
) {
  const int n_int = cum_p.size();
  NumericVector out(n);

  // Draw n uniforms via dqrng
  NumericVector u = dqrng::dqrunif(n);

  for (int i = 0; i < n; ++i) {
    const double ui = u[i];

    // Binary search for interval: smallest k such that cum_p[k] >= ui
    int lo = 0, hi = n_int - 1;
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      if (cum_p[mid] < ui) lo = mid + 1; else hi = mid;
    }
    const int k = lo;

    // Rescale ui to [0,1) within the selected interval
    const double p_lo = (k == 0) ? 0.0 : cum_p[k - 1];
    const double p_hi = cum_p[k];
    const double frac = (ui - p_lo) / (p_hi - p_lo);

    out[i] = a_time[k] + frac * (a_time[k + 1] - a_time[k]);
  }

  return out;
}

// ------------------------------------------------------------------ //
//  Internal: piecewise exponential sampler (survival / dropout times)
// ------------------------------------------------------------------ //
// Equivalent to R's rpiece_exp() for the piecewise case.
// Draws one Exp(1) variate per subject and maps it through the inverse
// cumulative hazard function in a single C++ loop, replacing
// findInterval(), pmax/pmin, and vector index operations.
//
// [[Rcpp::export]]
NumericVector rpiece_exp_cpp(
    int n,
    const NumericVector& hazard,   // hazard rates per interval
    const NumericVector& fin_time, // finite left endpoints (length = n_int)
    const NumericVector& cum_haz   // cumulative hazard at each left endpoint
) {
  const int n_int = hazard.size();
  NumericVector out(n);

  // Draw n Exp(1) variates via dqrng
  NumericVector target = dqrng::dqrexp(n, 1.0);

  for (int i = 0; i < n; ++i) {
    const double tgt = target[i];

    // Binary search: largest k such that cum_haz[k] <= tgt
    int lo = 0, hi = n_int - 1;
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      if (cum_haz[mid] <= tgt) lo = mid + 1; else hi = mid;
    }
    int k = lo;
    if (cum_haz[k] > tgt && k > 0) --k;

    out[i] = fin_time[k] + (tgt - cum_haz[k]) / hazard[k];
  }

  return out;
}

// ------------------------------------------------------------------ //
//  Internal: build sorted two-group output without rbind + order
// ------------------------------------------------------------------ //
// Interleaves control (group=1) and treatment (group=2) rows in sim order
// directly, avoiding rbind() followed by order(df$sim, df$group).
// Uses memcpy for contiguous blocks of rows (per simulation, per group),
// which is substantially faster than per-element copies.
//
// [[Rcpp::export]]
DataFrame interleave_groups(
    const DataFrame& df_ctrl,
    const DataFrame& df_trt,
    int nsim,
    int n_ctrl,
    int n_trt
) {
  const int total = nsim * (n_ctrl + n_trt);

  IntegerVector  sim_out(total);
  IntegerVector  grp_out(total);
  NumericVector  acc_out(total), srv_out(total), dro_out(total);
  NumericVector  tte_out(total), cal_out(total);
  IntegerVector  evt_out(total);

  // Extract column references once (cheap: just SEXP pointer assignment)
  IntegerVector sim_c  = df_ctrl["sim"];
  IntegerVector grp_c  = df_ctrl["group"];
  NumericVector acc_c  = df_ctrl["accrual_time"];
  NumericVector srv_c  = df_ctrl["surv_time"];
  NumericVector dro_c  = df_ctrl["dropout_time"];
  NumericVector tte_c  = df_ctrl["tte"];
  IntegerVector evt_c  = df_ctrl["event"];
  NumericVector cal_c  = df_ctrl["calendar_time"];

  IntegerVector sim_t  = df_trt["sim"];
  IntegerVector grp_t  = df_trt["group"];
  NumericVector acc_t  = df_trt["accrual_time"];
  NumericVector srv_t  = df_trt["surv_time"];
  NumericVector dro_t  = df_trt["dropout_time"];
  NumericVector tte_t  = df_trt["tte"];
  IntegerVector evt_t  = df_trt["event"];
  NumericVector cal_t  = df_trt["calendar_time"];

  // Raw pointers for memcpy
  int*    p_sim_out = INTEGER(sim_out);
  int*    p_grp_out = INTEGER(grp_out);
  double* p_acc_out = REAL(acc_out);
  double* p_srv_out = REAL(srv_out);
  double* p_dro_out = REAL(dro_out);
  double* p_tte_out = REAL(tte_out);
  int*    p_evt_out = INTEGER(evt_out);
  double* p_cal_out = REAL(cal_out);

  const int*    p_sim_c = INTEGER(sim_c);
  const int*    p_grp_c = INTEGER(grp_c);
  const double* p_acc_c = REAL(acc_c);
  const double* p_srv_c = REAL(srv_c);
  const double* p_dro_c = REAL(dro_c);
  const double* p_tte_c = REAL(tte_c);
  const int*    p_evt_c = INTEGER(evt_c);
  const double* p_cal_c = REAL(cal_c);

  const int*    p_sim_t = INTEGER(sim_t);
  const int*    p_grp_t = INTEGER(grp_t);
  const double* p_acc_t = REAL(acc_t);
  const double* p_srv_t = REAL(srv_t);
  const double* p_dro_t = REAL(dro_t);
  const double* p_tte_t = REAL(tte_t);
  const int*    p_evt_t = INTEGER(evt_t);
  const double* p_cal_t = REAL(cal_t);

  int out_idx = 0;
  for (int s = 0; s < nsim; ++s) {
    // Control block (n_ctrl rows)
    const int c_base = s * n_ctrl;
    std::memcpy(p_sim_out + out_idx, p_sim_c + c_base, n_ctrl * sizeof(int));
    std::memcpy(p_grp_out + out_idx, p_grp_c + c_base, n_ctrl * sizeof(int));
    std::memcpy(p_acc_out + out_idx, p_acc_c + c_base, n_ctrl * sizeof(double));
    std::memcpy(p_srv_out + out_idx, p_srv_c + c_base, n_ctrl * sizeof(double));
    std::memcpy(p_dro_out + out_idx, p_dro_c + c_base, n_ctrl * sizeof(double));
    std::memcpy(p_tte_out + out_idx, p_tte_c + c_base, n_ctrl * sizeof(double));
    std::memcpy(p_evt_out + out_idx, p_evt_c + c_base, n_ctrl * sizeof(int));
    std::memcpy(p_cal_out + out_idx, p_cal_c + c_base, n_ctrl * sizeof(double));
    out_idx += n_ctrl;

    // Treatment block (n_trt rows)
    const int t_base = s * n_trt;
    std::memcpy(p_sim_out + out_idx, p_sim_t + t_base, n_trt * sizeof(int));
    std::memcpy(p_grp_out + out_idx, p_grp_t + t_base, n_trt * sizeof(int));
    std::memcpy(p_acc_out + out_idx, p_acc_t + t_base, n_trt * sizeof(double));
    std::memcpy(p_srv_out + out_idx, p_srv_t + t_base, n_trt * sizeof(double));
    std::memcpy(p_dro_out + out_idx, p_dro_t + t_base, n_trt * sizeof(double));
    std::memcpy(p_tte_out + out_idx, p_tte_t + t_base, n_trt * sizeof(double));
    std::memcpy(p_evt_out + out_idx, p_evt_t + t_base, n_trt * sizeof(int));
    std::memcpy(p_cal_out + out_idx, p_cal_t + t_base, n_trt * sizeof(double));
    out_idx += n_trt;
  }

  return DataFrame::create(
    Named("sim")           = sim_out,
    Named("group")         = grp_out,
    Named("accrual_time")  = acc_out,
    Named("surv_time")     = srv_out,
    Named("dropout_time")  = dro_out,
    Named("tte")           = tte_out,
    Named("event")         = evt_out,
    Named("calendar_time") = cal_out
  );
}
