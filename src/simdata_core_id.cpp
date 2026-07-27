// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(dqrng)]]
#include <Rcpp.h>
#include <dqrng.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ------------------------------------------------------------------ //
//  Illness-death (three-state) simulator for two correlated
//  time-to-event endpoints with optional treatment switching
// ------------------------------------------------------------------ //
// This kernel is a separate entry point from simdata_core_full and shares none
// of its state. It generates two correlated endpoints from the standard
// illness-death structure with three transition hazards:
//
//   h_01 : state 0 -> 1, the non-terminal (intermediate) event
//   h_02 : state 0 -> 2, the terminal event without an intermediate event
//   h_12 : state 1 -> 2, the terminal event after an intermediate event
//          (no-switch or, for switchers, a separate post-switch hazard)
//
// Two latent times are drawn per subject: T01 from h_01 and T02 from h_02. The
// first endpoint is the state-0 exit time min(T01, T02); an intermediate event
// occurs when T01 < T02. A direct terminal event (T02 <= T01) gives terminal
// time T02. After an intermediate event the post-event survival OS' is drawn
// from h_12 on a clock-reset (semi-Markov) basis and the terminal time is
// T01 + OS'. With h_12 equal to h_02 and switch probability zero this reduces
// exactly to the Fleischer (2009) maximal-independence model (their Theorem 1).
//
// In oncology the first endpoint is progression-free survival, the terminal
// endpoint is overall survival, and the intermediate event is progression; the
// same structure applies to other settings (for example a non-fatal event and
// all-cause death in a cardiovascular outcome trial). Subgroups are not
// supported here; the no-subgroup pipeline is used.

// Per-group dqrng consumption order (drawn contiguously over all simulations):
//   1. accrual (one dqrunif per subject, within-interval position)
//   2. T01     (one dqrexp per subject, intermediate-event latent time)
//   3. T02     (one dqrexp per subject, direct terminal latent time)
//   4. switch  (one dqrunif per subject, switch vs no-switch given an event)
//   5. OS'     (one dqrexp per subject, post-event survival)
//   6. dropout (one dqrexp per subject, only when has_dropout is true)
// In two-group mode the whole control block is drawn before the whole treatment
// block; rows are interleaved into (sim, group) order afterwards (no RNG cost).

// Deterministic piecewise accrual: each interval receives a fixed count
// (supplied by the caller, summing to the group size), repeated every
// simulation, with uniform within-interval entry positions. One dqrunif per
// subject in subject order.
static void draw_accrual_det(double* out, int nsim, int n,
                             const std::vector<double>& a_time,
                             const std::vector<int>& counts) {
  const int n_int   = (int) counts.size();
  const int total_n = nsim * n;
  NumericVector u = dqrng::dqrunif(total_n);
  for (int s = 0; s < nsim; ++s) {
    const int base_s = s * n;
    int j = 0;
    for (int i = 0; i < n_int; ++i) {
      const double lo    = a_time[i];
      const double width = a_time[i + 1] - a_time[i];
      const int    ci    = counts[i];
      for (int r = 0; r < ci; ++r, ++j) {
        const int gi = base_s + j;
        out[gi] = lo + u[gi] * width;
      }
    }
  }
}

// Invert a unit-exponential target through a piecewise-exponential cumulative
// hazard. For a single hazard this is target / hazard. Mirrors the inversion
// used by draw_exp in simdata_core_full.
static inline double inv_cum_hazard(double tgt,
                                    const std::vector<double>& hazard,
                                    const std::vector<double>& fin_time,
                                    const std::vector<double>& cum_haz) {
  const int n_int = (int) hazard.size();
  if (n_int == 1) return tgt / hazard[0];
  int lo = 0, hi = n_int - 1;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (cum_haz[mid] <= tgt) lo = mid + 1; else hi = mid;
  }
  int k = lo;
  if (cum_haz[k] > tgt && k > 0) --k;
  return fin_time[k] + (tgt - cum_haz[k]) / hazard[k];
}

// Fill one group's block of every output column starting at row 'base'.
static void simulate_group_id_into(
    int base, int nsim, int n, int group_id,
    const std::vector<double>& a_time, const std::vector<int>& acc_counts,
    const std::vector<double>& h01_haz, const std::vector<double>& h01_fin, const std::vector<double>& h01_cum,
    const std::vector<double>& h02_haz, const std::vector<double>& h02_fin, const std::vector<double>& h02_cum,
    const std::vector<double>& h12_haz, const std::vector<double>& h12_fin, const std::vector<double>& h12_cum,
    const std::vector<double>& h12s_haz, const std::vector<double>& h12s_fin, const std::vector<double>& h12s_cum,
    double switch_prop,
    bool has_dropout,
    const std::vector<double>& d_haz, const std::vector<double>& d_fin, const std::vector<double>& d_cum,
    int* sim_col, int* grp_col,
    double* acc_col, double* e1_surv_col, double* e2_surv_col, double* dro_col,
    double* e1_tte_col, int* e1_evt_col,
    double* e2_tte_col, int* e2_evt_col,
    double* e1_cal_col, double* e2_cal_col,
    int* inter_col, int* sw_col, double* sw_time_col
) {
  const int total_n = nsim * n;

  // sim and group columns.
  for (int s = 0; s < nsim; ++s)
    for (int k = 0; k < n; ++k)
      sim_col[base + s * n + k] = s + 1;
  for (int i = 0; i < total_n; ++i) grp_col[base + i] = group_id;

  // 1. accrual.
  draw_accrual_det(acc_col + base, nsim, n, a_time, acc_counts);

  // 2-5. transition-process draws (whole block, fixed order).
  NumericVector e_01  = dqrng::dqrexp(total_n, 1.0);
  NumericVector e_02  = dqrng::dqrexp(total_n, 1.0);
  NumericVector u_sw  = dqrng::dqrunif(total_n);
  NumericVector e_osp = dqrng::dqrexp(total_n, 1.0);

  // 6. dropout (shared by both endpoints) or +Inf.
  if (has_dropout) {
    NumericVector e_dro = dqrng::dqrexp(total_n, 1.0);
    for (int i = 0; i < total_n; ++i)
      dro_col[base + i] = inv_cum_hazard(e_dro[i], d_haz, d_fin, d_cum);
  } else {
    for (int i = 0; i < total_n; ++i) dro_col[base + i] = R_PosInf;
  }

  for (int i = 0; i < total_n; ++i) {
    const double t01 = inv_cum_hazard(e_01[i], h01_haz, h01_fin, h01_cum);
    const double t02 = inv_cum_hazard(e_02[i], h02_haz, h02_fin, h02_cum);

    double e1_t, e2_t;
    int    intermediate = 0, switched = 0;
    double sw_time = NA_REAL;

    if (t02 <= t01) {
      // Direct terminal event before any intermediate event.
      e1_t = t02;
      e2_t = t02;
    } else {
      // Intermediate event first, then post-event survival (clock-reset).
      intermediate = 1;
      e1_t = t01;
      const bool sw = (u_sw[i] < switch_prop);
      double osp;
      if (sw) {
        switched = 1;
        sw_time  = t01;
        osp = inv_cum_hazard(e_osp[i], h12s_haz, h12s_fin, h12s_cum);
      } else {
        osp = inv_cum_hazard(e_osp[i], h12_haz, h12_fin, h12_cum);
      }
      e2_t = t01 + osp;
    }

    const double dr = dro_col[base + i];
    const double e1_tt = (e1_t <= dr) ? e1_t : dr;
    const double e2_tt = (e2_t <= dr) ? e2_t : dr;

    e1_surv_col[base + i] = e1_t;
    e2_surv_col[base + i] = e2_t;
    e1_tte_col[base + i]  = e1_tt;
    e1_evt_col[base + i]  = (e1_t <= dr) ? 1 : 0;
    e2_tte_col[base + i]  = e2_tt;
    e2_evt_col[base + i]  = (e2_t <= dr) ? 1 : 0;
    e1_cal_col[base + i]  = acc_col[base + i] + e1_tt;
    e2_cal_col[base + i]  = acc_col[base + i] + e2_tt;
    inter_col[base + i]   = intermediate;
    sw_col[base + i]      = switched;
    sw_time_col[base + i] = sw_time;
  }
}

// ------------------------------------------------------------------ //
//  Whole-pipeline illness-death simulator
// ------------------------------------------------------------------ //
// [[Rcpp::export]]
DataFrame simdata_core_id(
    int nsim,
    const IntegerVector& n_grp,          // length 1 or 2
    const NumericVector& a_time,
    const IntegerVector& acc_counts_c,
    const IntegerVector& acc_counts_t,
    const NumericVector& h01_haz_c, const NumericVector& h01_fin_c, const NumericVector& h01_cum_c,
    const NumericVector& h01_haz_t, const NumericVector& h01_fin_t, const NumericVector& h01_cum_t,
    const NumericVector& h02_haz_c, const NumericVector& h02_fin_c, const NumericVector& h02_cum_c,
    const NumericVector& h02_haz_t, const NumericVector& h02_fin_t, const NumericVector& h02_cum_t,
    const NumericVector& h12_haz_c, const NumericVector& h12_fin_c, const NumericVector& h12_cum_c,
    const NumericVector& h12_haz_t, const NumericVector& h12_fin_t, const NumericVector& h12_cum_t,
    const NumericVector& h12s_haz_c, const NumericVector& h12s_fin_c, const NumericVector& h12s_cum_c,
    const NumericVector& h12s_haz_t, const NumericVector& h12s_fin_t, const NumericVector& h12s_cum_t,
    double switch_prop_c, double switch_prop_t,
    bool has_dropout,
    const NumericVector& d_haz_c, const NumericVector& d_fin_c, const NumericVector& d_cum_c,
    const NumericVector& d_haz_t, const NumericVector& d_fin_t, const NumericVector& d_cum_t
) {
  const int n_groups = n_grp.size();
  const std::vector<double> at(a_time.begin(), a_time.end());
  const std::vector<int> acc_c(acc_counts_c.begin(), acc_counts_c.end());
  const std::vector<int> acc_t(acc_counts_t.begin(), acc_counts_t.end());

  auto v = [](const NumericVector& x) {
    return std::vector<double>(x.begin(), x.end());
  };
  const auto h01h_c = v(h01_haz_c), h01f_c = v(h01_fin_c), h01c_c = v(h01_cum_c);
  const auto h01h_t = v(h01_haz_t), h01f_t = v(h01_fin_t), h01c_t = v(h01_cum_t);
  const auto h02h_c = v(h02_haz_c), h02f_c = v(h02_fin_c), h02c_c = v(h02_cum_c);
  const auto h02h_t = v(h02_haz_t), h02f_t = v(h02_fin_t), h02c_t = v(h02_cum_t);
  const auto h12h_c = v(h12_haz_c), h12f_c = v(h12_fin_c), h12c_c = v(h12_cum_c);
  const auto h12h_t = v(h12_haz_t), h12f_t = v(h12_fin_t), h12c_t = v(h12_cum_t);
  const auto h12sh_c = v(h12s_haz_c), h12sf_c = v(h12s_fin_c), h12sc_c = v(h12s_cum_c);
  const auto h12sh_t = v(h12s_haz_t), h12sf_t = v(h12s_fin_t), h12sc_t = v(h12s_cum_t);
  const auto dh_c = v(d_haz_c), df_c = v(d_fin_c), dc_c = v(d_cum_c);
  const auto dh_t = v(d_haz_t), df_t = v(d_fin_t), dc_t = v(d_cum_t);

  auto assemble = [](IntegerVector sim_col, IntegerVector grp_col,
                     NumericVector acc_col, NumericVector e1_surv_col,
                     NumericVector e2_surv_col, NumericVector dro_col,
                     NumericVector e1_tte_col, IntegerVector e1_evt_col,
                     NumericVector e2_tte_col, IntegerVector e2_evt_col,
                     NumericVector e1_cal_col, NumericVector e2_cal_col,
                     IntegerVector inter_col, IntegerVector sw_col,
                     NumericVector sw_time_col) {
    List cols;
    cols["sim"]              = sim_col;
    cols["group"]            = grp_col;
    cols["accrual_time"]     = acc_col;
    cols["e1_surv_time"]     = e1_surv_col;
    cols["e2_surv_time"]     = e2_surv_col;
    cols["dropout_time"]     = dro_col;
    cols["e1_tte"]           = e1_tte_col;
    cols["e1_event"]         = e1_evt_col;
    cols["e2_tte"]           = e2_tte_col;
    cols["e2_event"]         = e2_evt_col;
    cols["e1_calendar_time"] = e1_cal_col;
    cols["e2_calendar_time"] = e2_cal_col;
    cols["intermediate"]     = inter_col;
    cols["switched"]         = sw_col;
    cols["switch_time"]      = sw_time_col;
    return DataFrame(cols);
  };

  if (n_groups == 1) {
    const int n0 = n_grp[0];
    const int total = nsim * n0;
    IntegerVector sim_col(total), grp_col(total), e1_evt(total), e2_evt(total),
                  inter(total), sw(total);
    NumericVector acc(total), e1_surv(total), e2_surv(total), dro(total),
                  e1_tte(total), e2_tte(total), e1_cal(total), e2_cal(total),
                  sw_time(total);

    simulate_group_id_into(
      0, nsim, n0, 1, at, acc_c,
      h01h_c, h01f_c, h01c_c, h02h_c, h02f_c, h02c_c,
      h12h_c, h12f_c, h12c_c, h12sh_c, h12sf_c, h12sc_c,
      switch_prop_c, has_dropout, dh_c, df_c, dc_c,
      INTEGER(sim_col), INTEGER(grp_col),
      REAL(acc), REAL(e1_surv), REAL(e2_surv), REAL(dro),
      REAL(e1_tte), INTEGER(e1_evt), REAL(e2_tte), INTEGER(e2_evt),
      REAL(e1_cal), REAL(e2_cal),
      INTEGER(inter), INTEGER(sw), REAL(sw_time));

    return assemble(sim_col, grp_col, acc, e1_surv, e2_surv, dro,
                    e1_tte, e1_evt, e2_tte, e2_evt, e1_cal, e2_cal,
                    inter, sw, sw_time);
  }

  // Two groups: draw each group contiguously (preserving the dqrng stream),
  // then interleave into (sim, group) order.
  const int nc = n_grp[0], nt = n_grp[1];
  const int per_sim = nc + nt;
  const int total = nsim * per_sim;
  const int totc = nsim * nc, tott = nsim * nt;

  IntegerVector simC(totc), grpC(totc), e1eC(totc), e2eC(totc), inC(totc), swC(totc);
  NumericVector acC(totc), e1sC(totc), e2sC(totc), drC(totc), e1tC(totc), e2tC(totc),
                e1cC(totc), e2cC(totc), stC(totc);
  IntegerVector simT(tott), grpT(tott), e1eT(tott), e2eT(tott), inT(tott), swT(tott);
  NumericVector acT(tott), e1sT(tott), e2sT(tott), drT(tott), e1tT(tott), e2tT(tott),
                e1cT(tott), e2cT(tott), stT(tott);

  simulate_group_id_into(
    0, nsim, nc, 1, at, acc_c,
    h01h_c, h01f_c, h01c_c, h02h_c, h02f_c, h02c_c,
    h12h_c, h12f_c, h12c_c, h12sh_c, h12sf_c, h12sc_c,
    switch_prop_c, has_dropout, dh_c, df_c, dc_c,
    INTEGER(simC), INTEGER(grpC),
    REAL(acC), REAL(e1sC), REAL(e2sC), REAL(drC),
    REAL(e1tC), INTEGER(e1eC), REAL(e2tC), INTEGER(e2eC),
    REAL(e1cC), REAL(e2cC), INTEGER(inC), INTEGER(swC), REAL(stC));

  simulate_group_id_into(
    0, nsim, nt, 2, at, acc_t,
    h01h_t, h01f_t, h01c_t, h02h_t, h02f_t, h02c_t,
    h12h_t, h12f_t, h12c_t, h12sh_t, h12sf_t, h12sc_t,
    switch_prop_t, has_dropout, dh_t, df_t, dc_t,
    INTEGER(simT), INTEGER(grpT),
    REAL(acT), REAL(e1sT), REAL(e2sT), REAL(drT),
    REAL(e1tT), INTEGER(e1eT), REAL(e2tT), INTEGER(e2eT),
    REAL(e1cT), REAL(e2cT), INTEGER(inT), INTEGER(swT), REAL(stT));

  IntegerVector sim_col(total), grp_col(total), e1_evt(total), e2_evt(total),
                inter(total), sw(total);
  NumericVector acc(total), e1_surv(total), e2_surv(total), dro(total),
                e1_tte(total), e2_tte(total), e1_cal(total), e2_cal(total),
                sw_time(total);

  int* p_sim = INTEGER(sim_col); int* p_grp = INTEGER(grp_col);
  int* p_e1e = INTEGER(e1_evt); int* p_e2e = INTEGER(e2_evt);
  int* p_in = INTEGER(inter); int* p_sw = INTEGER(sw);
  double* p_ac = REAL(acc); double* p_e1s = REAL(e1_surv); double* p_e2s = REAL(e2_surv);
  double* p_dr = REAL(dro); double* p_e1t = REAL(e1_tte); double* p_e2t = REAL(e2_tte);
  double* p_e1c = REAL(e1_cal); double* p_e2c = REAL(e2_cal); double* p_st = REAL(sw_time);

  const int* sC = INTEGER(simC); const int* gC = INTEGER(grpC);
  const int* e1eCp = INTEGER(e1eC); const int* e2eCp = INTEGER(e2eC);
  const int* inCp = INTEGER(inC); const int* swCp = INTEGER(swC);
  const double* acCp = REAL(acC); const double* e1sCp = REAL(e1sC); const double* e2sCp = REAL(e2sC);
  const double* drCp = REAL(drC); const double* e1tCp = REAL(e1tC); const double* e2tCp = REAL(e2tC);
  const double* e1cCp = REAL(e1cC); const double* e2cCp = REAL(e2cC); const double* stCp = REAL(stC);
  const int* sT = INTEGER(simT); const int* gT = INTEGER(grpT);
  const int* e1eTp = INTEGER(e1eT); const int* e2eTp = INTEGER(e2eT);
  const int* inTp = INTEGER(inT); const int* swTp = INTEGER(swT);
  const double* acTp = REAL(acT); const double* e1sTp = REAL(e1sT); const double* e2sTp = REAL(e2sT);
  const double* drTp = REAL(drT); const double* e1tTp = REAL(e1tT); const double* e2tTp = REAL(e2tT);
  const double* e1cTp = REAL(e1cT); const double* e2cTp = REAL(e2cT); const double* stTp = REAL(stT);

  int out = 0;
  for (int s = 0; s < nsim; ++s) {
    const int cb = s * nc;
    for (int k = 0; k < nc; ++k, ++out) {
      p_sim[out] = sC[cb + k]; p_grp[out] = gC[cb + k];
      p_e1e[out] = e1eCp[cb + k]; p_e2e[out] = e2eCp[cb + k];
      p_in[out] = inCp[cb + k]; p_sw[out] = swCp[cb + k];
      p_ac[out] = acCp[cb + k]; p_e1s[out] = e1sCp[cb + k]; p_e2s[out] = e2sCp[cb + k];
      p_dr[out] = drCp[cb + k]; p_e1t[out] = e1tCp[cb + k]; p_e2t[out] = e2tCp[cb + k];
      p_e1c[out] = e1cCp[cb + k]; p_e2c[out] = e2cCp[cb + k]; p_st[out] = stCp[cb + k];
    }
    const int tb = s * nt;
    for (int k = 0; k < nt; ++k, ++out) {
      p_sim[out] = sT[tb + k]; p_grp[out] = gT[tb + k];
      p_e1e[out] = e1eTp[tb + k]; p_e2e[out] = e2eTp[tb + k];
      p_in[out] = inTp[tb + k]; p_sw[out] = swTp[tb + k];
      p_ac[out] = acTp[tb + k]; p_e1s[out] = e1sTp[tb + k]; p_e2s[out] = e2sTp[tb + k];
      p_dr[out] = drTp[tb + k]; p_e1t[out] = e1tTp[tb + k]; p_e2t[out] = e2tTp[tb + k];
      p_e1c[out] = e1cTp[tb + k]; p_e2c[out] = e2cTp[tb + k]; p_st[out] = stTp[tb + k];
    }
  }

  return assemble(sim_col, grp_col, acc, e1_surv, e2_surv, dro,
                  e1_tte, e1_evt, e2_tte, e2_evt, e1_cal, e2_cal,
                  inter, sw, sw_time);
}
