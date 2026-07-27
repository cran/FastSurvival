// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
using namespace Rcpp;

// Forward declarations of the pointer-based statistic cores. Each is an
// ordinary C++ function with external linkage defined in its own source file
// (the exported Rcpp wrappers in those files now delegate to these). Resolving
// them here uses the same cross-translation-unit linkage RcppExports.cpp
// relies on. Signatures must match the definitions exactly.
void logrank_core_impl(const double*, const int*, const int*, int, double*);
void weighted_logrank_core_impl(const double*, const int*, const int*, int,
                                int, double, double, double, double*);
void stratified_logrank_core_impl(const double*, const int*, const int*,
                                  const int*, int, double*);
void stratified_weighted_logrank_core_impl(const double*, const int*,
                                           const int*, const int*, int,
                                           int, double, double, double,
                                           double*);
struct PiheEvSummary { double n1k; double n0k; double OTk; double Ok; };
void pihe_core_impl(const double*, const int*, const int*, int,
                    std::vector<PiheEvSummary>&, double*);
struct RmstEvSummary { double cum_area; double g; };
void rmst_core_impl(const double*, const double*, int, double,
                    std::vector<RmstEvSummary>&, double*);
void km_core_impl(const double*, const double*, int, double, double*);
double combo_logrank_core_impl(const double*, const int*, const int*, int,
                               const double*, const double*, int,
                               double*, double*, std::vector<double>&);
void ahsw_core_impl(const double*, const int*, int, double,
                    std::vector<double>&, std::vector<double>&,
                    std::vector<double>&, double*);
void milestone_core_impl(const double*, const int*, const int*, int, double,
                         double*);
void rmw_core_impl(const double*, const int*, const int*, int, double, double*);
void ahr_core_two_impl(const double*, const int*, int, const double*,
                       const int*, int, double, double*);
struct WmstEvSummary {
  double e_time; double at_risk; int n_event; double surv_after; double tail;
};
void wmst_core_impl(const double*, const int*, int, double, double,
                    std::vector<WmstEvSummary>&, double*);
struct WkmScratch {
  std::vector<double> width; std::vector<double> wt;
  std::vector<double> s1; std::vector<double> s2;
  std::vector<double> sp; std::vector<double> a_seq;
};
void wkm_core_impl(const double*, const int*, const int*, int, int,
                   WkmScratch&, double*);
struct MedsurvScratch {
  std::vector<double> event_time; std::vector<int> event_count;
  std::vector<double> at_risk; std::vector<double> hazard_incr;
};
void medsurv_core_impl(const double*, const int*, int, double,
                       MedsurvScratch&, double*);

// Silverman type default bandwidth for the median kernel hazard, computed from
// the event times of a single group. Reproduces the default bw used by
// medsurv_fast when bw is not supplied: 1.06 sd(evt) m^(-1/5) when the event
// times have positive spread, falling back to range / max(m, 1)^(1/3), or 1.
static double silverman_bw(const double* t, const int* e, int n) {
  double sum = 0.0;
  int m = 0;
  for (int k = 0; k < n; ++k) {
    if (e[k] == 1) { sum += t[k]; ++m; }
  }
  if (m >= 2) {
    double mean = sum / (double) m;
    double ss = 0.0;
    for (int k = 0; k < n; ++k) {
      if (e[k] == 1) { double d = t[k] - mean; ss += d * d; }
    }
    double sd = std::sqrt(ss / (double) (m - 1));
    if (sd > 0.0) {
      double b = 1.06 * sd * std::pow((double) m, -1.0 / 5.0);
      if (std::isfinite(b) && b > 0.0) return b;
    }
  }
  double tmin = R_PosInf, tmax = R_NegInf;
  for (int k = 0; k < n; ++k) {
    if (e[k] == 1) {
      if (t[k] < tmin) tmin = t[k];
      if (t[k] > tmax) tmax = t[k];
    }
  }
  double span = (m >= 1) ? (tmax - tmin) : 0.0;
  double mm = (m >= 1) ? (double) m : 1.0;
  if (span > 0.0) return span / std::pow(mm, 1.0 / 3.0);
  return 1.0;
}

// Stable LSD radix sort of indices [0, m) by a non-negative double key. The
// observed time at an administrative cut is always non-negative, so the IEEE
// 754 bit pattern of a non-negative double is monotone in the value and can be
// sorted as a 64-bit unsigned integer. The sort is stable, so tied times keep
// their input order, which matches the std::sort comparator the kernel used
// before and the atomic tied-block handling in every statistic core. idx is
// filled with the sort order; key, scratch, and cnt are reusable buffers.
static void radix_sort_idx(const std::vector<double>& key, int m,
                           std::vector<int>& idx, std::vector<int>& scratch,
                           std::vector<uint64_t>& bits) {
  if ((int) idx.size() < m) idx.resize(m);
  if ((int) scratch.size() < m) scratch.resize(m);
  if ((int) bits.size() < m) bits.resize(m);

  for (int k = 0; k < m; ++k) {
    double v = key[k];
    uint64_t u;
    std::memcpy(&u, &v, sizeof(u));   // non-negative double: bits are monotone
    bits[k] = u;
    idx[k] = k;
  }

  int* in  = idx.data();
  int* out = scratch.data();
  int cnt[256];

  for (int shift = 0; shift < 64; shift += 8) {
    std::memset(cnt, 0, sizeof(cnt));
    for (int k = 0; k < m; ++k) {
      ++cnt[(bits[in[k]] >> shift) & 0xFF];
    }
    int sum = 0;
    for (int b = 0; b < 256; ++b) { int c = cnt[b]; cnt[b] = sum; sum += c; }
    for (int k = 0; k < m; ++k) {
      const int v = in[k];
      out[cnt[(bits[v] >> shift) & 0xFF]++] = v;
    }
    int* t = in; in = out; out = t;
  }

  // After 8 passes (even count) the sorted order is back in idx.data(); guard
  // anyway in case the buffers were swapped an odd number of times.
  if (in != idx.data()) {
    std::memcpy(idx.data(), in, sizeof(int) * (size_t) m);
  }
}

// Fused analysis loop. Walks simulations, looks, and populations in C++,
// applies the administrative cut inline, and calls the pointer-based statistic
// cores directly on reusable buffers, so no Rcpp object is allocated per
// (sim, look, population) cell. Returns the raw core outputs; all conversion
// into z-scores, confidence intervals, and p-values (and the max-combo
// multivariate-normal p-value) is done in the R wrapper.
//
// The data must be sorted by simulation, with sim_ptr giving 0-based row
// offsets so that simulation si occupies rows [sim_ptr[si], sim_ptr[si + 1]).
// Output cells are ordered (sim, look, population) to match analysis_fast.
//
// [[Rcpp::export]]
List analysis_loop_core(
    const IntegerVector& sim_ptr,
    const NumericVector& accrual,
    const NumericVector& tte,
    const IntegerVector& event,
    const IntegerVector& j,
    int look_type,                  // 0 = event-driven, 1 = calendar time
    const NumericVector& look_values,
    const IntegerVector& pop_col,   // -1 = overall, else 0-based column in sub_mat
    const IntegerVector& pop_level,
    const IntegerMatrix& sub_mat,   // N x n_subcols (subgroup labels)
    const IntegerVector& strata,    // length N or length 0
    bool use_strata,
    bool do_logrank, bool do_coxph, bool do_rmst, bool do_km,
    bool do_maxcombo, bool do_ahsw,
    bool do_milestone, bool do_rmw, bool do_ahr,
    bool do_medsurv, bool do_wkm, bool do_wmst,
    int weight_scheme,              // -1 plain, 0 FH, 1 mwlrt, 2 gehan, 3 tarone-ware
    double rho, double gamma, double t_star,
    const NumericVector& mc_rho, const NumericVector& mc_gamma,
    double tau, double t_eval, double s_star,
    double wmst_tau1, double wmst_tau2, int wkm_weight, double medsurv_bw
) {
  const int nsim    = sim_ptr.size() - 1;
  const int n_looks = look_values.size();
  const int n_pop   = pop_col.size();
  const int ncells  = nsim * n_looks * n_pop;
  const int nw      = mc_rho.size();

  NumericVector cutoff_out(ncells, NA_REAL);
  IntegerVector reached_out(ncells);
  IntegerVector nen_out(ncells);
  IntegerVector nev_out(ncells);
  IntegerVector nd_out(ncells);

  NumericMatrix lr_mat, cox_mat, rmst_mat, km_mat, ahsw_mat, mc_Uout, mc_Vout;
  NumericMatrix ms_mat, rmw_mat, ahr_mat;
  NumericMatrix medsurv_mat, wkm_mat, wmst_mat;
  if (do_logrank) { lr_mat   = NumericMatrix(ncells, 2); std::fill(lr_mat.begin(),   lr_mat.end(),   NA_REAL); }
  if (do_coxph)   { cox_mat  = NumericMatrix(ncells, 4); std::fill(cox_mat.begin(),  cox_mat.end(),  NA_REAL); }
  if (do_rmst)    { rmst_mat = NumericMatrix(ncells, 4); std::fill(rmst_mat.begin(), rmst_mat.end(), NA_REAL); }
  if (do_km)      { km_mat   = NumericMatrix(ncells, 2); std::fill(km_mat.begin(),   km_mat.end(),   NA_REAL); }
  if (do_ahsw)    { ahsw_mat = NumericMatrix(ncells, 8); std::fill(ahsw_mat.begin(), ahsw_mat.end(), NA_REAL); }
  if (do_milestone) { ms_mat  = NumericMatrix(ncells, 4); std::fill(ms_mat.begin(),  ms_mat.end(),  NA_REAL); }
  if (do_rmw)       { rmw_mat = NumericMatrix(ncells, 5); std::fill(rmw_mat.begin(), rmw_mat.end(), NA_REAL); }
  if (do_ahr)       { ahr_mat = NumericMatrix(ncells, 5); std::fill(ahr_mat.begin(), ahr_mat.end(), NA_REAL); }
  if (do_medsurv) { medsurv_mat = NumericMatrix(ncells, 10); std::fill(medsurv_mat.begin(), medsurv_mat.end(), NA_REAL); }
  if (do_wkm)     { wkm_mat  = NumericMatrix(ncells, 4); std::fill(wkm_mat.begin(),  wkm_mat.end(),  NA_REAL); }
  if (do_wmst)    { wmst_mat = NumericMatrix(ncells, 4); std::fill(wmst_mat.begin(), wmst_mat.end(), NA_REAL); }
  if (do_maxcombo) {
    mc_Uout = NumericMatrix(ncells, nw);      std::fill(mc_Uout.begin(), mc_Uout.end(), NA_REAL);
    mc_Vout = NumericMatrix(ncells, nw * nw); std::fill(mc_Vout.begin(), mc_Vout.end(), NA_REAL);
  }

  // Largest simulation block, used to size the reusable buffers once.
  int max_block = 0;
  for (int si = 0; si < nsim; ++si) {
    const int sz = sim_ptr[si + 1] - sim_ptr[si];
    if (sz > max_block) max_block = sz;
  }

  // ---- Reusable buffers, allocated once outside the loops ----------------
  std::vector<double> cal_ev;            cal_ev.reserve(max_block);
  std::vector<double> cal_ev_tmp;        cal_ev_tmp.reserve(max_block);
  std::vector<int>    orig_cut(max_block);
  std::vector<double> t_cut(max_block);
  std::vector<int>    e_cut(max_block);
  std::vector<int>    j_cut(max_block);
  std::vector<int>    drop_cut(max_block);
  std::vector<int>    ord(max_block);
  std::vector<int>    radix_scratch(max_block);
  std::vector<uint64_t> radix_bits(max_block);

  // Population subset (overall reuses ord directly; subgroup fills sel).
  std::vector<int> sel(max_block);

  // Pooled population vectors (time-sorted), pointer-fed to the cores.
  std::vector<double> t_sel(max_block);
  std::vector<int>    ei_sel(max_block);
  std::vector<int>    j_sel(max_block);

  // Group-split buffers reused by rmst / km / ahsw / wmst / medsurv.
  std::vector<double> t0d(max_block), t1d(max_block);
  std::vector<int>    e0i(max_block), e1i(max_block);
  std::vector<double> e0d(max_block), e1d(max_block);

  // Stratified reorder buffers.
  std::vector<int>    st_sel(max_block), so(max_block);
  std::vector<double> ts(max_block);
  std::vector<int>    es(max_block), js(max_block), ss(max_block);

  // Core scratch reused across cells.
  std::vector<PiheEvSummary> pihe_ev;
  std::vector<RmstEvSummary> rmst_ev;
  std::vector<double> ahsw_dH, ahsw_Rrun, ahsw_Gfrac;
  std::vector<double> mc_w(nw > 0 ? nw : 1);
  std::vector<double> mc_ubuf((size_t)(nw > 0 ? nw : 1));
  std::vector<double> mc_vbuf((size_t)(nw > 0 ? nw : 1) * (nw > 0 ? nw : 1));
  std::vector<WmstEvSummary> wmst_ev;
  WkmScratch wkm_sc;
  MedsurvScratch medsurv_sc;

  int pos = 0;
  for (int si = 0; si < nsim; ++si) {
    const int g0 = sim_ptr[si];
    const int g1 = sim_ptr[si + 1];

    // Calendar event times for event-driven cutoffs (whole simulation).
    if (look_type == 0) {
      cal_ev.clear();
      for (int g = g0; g < g1; ++g) {
        if (event[g] == 1) cal_ev.push_back(accrual[g] + tte[g]);
      }
    }

    for (int l = 0; l < n_looks; ++l) {
      double cut_cutoff;
      bool   reached;
      if (look_type == 0) {
        const int target = (int) look_values[l];
        const int n_ev_sim = (int) cal_ev.size();
        if (target < 1 || target > n_ev_sim) {
          cut_cutoff = R_PosInf;
          reached    = false;
        } else {
          cal_ev_tmp.assign(cal_ev.begin(), cal_ev.end());
          std::nth_element(cal_ev_tmp.begin(),
                           cal_ev_tmp.begin() + (target - 1),
                           cal_ev_tmp.end());
          cut_cutoff = cal_ev_tmp[target - 1];
          reached    = true;
        }
      } else {
        cut_cutoff = look_values[l];
        reached    = true;
      }

      // Administrative cut for this (sim, look): enrolled subjects, censored
      // at the cutoff. orig_cut holds the within-sim 0-based index so subgroup
      // and stratum labels can be realigned.
      int m = 0;
      for (int g = g0; g < g1; ++g) {
        const double a = accrual[g];
        if (a > cut_cutoff) continue;
        const double full   = tte[g];
        const bool   before = (a + full <= cut_cutoff);
        orig_cut[m] = g - g0;
        t_cut[m]    = before ? full : (cut_cutoff - a);
        e_cut[m]    = before ? event[g] : 0;
        j_cut[m]    = j[g];
        // A dropout is an enrolled subject whose natural censoring (event == 0)
        // occurred on or before the cutoff; subjects not yet resolved at the
        // cutoff are administrative (pipeline), not dropouts.
        drop_cut[m] = (before && event[g] == 0) ? 1 : 0;
        ++m;
      }

      // Stable radix sort of the cut by observed time (non-negative keys).
      radix_sort_idx(t_cut, m, ord, radix_scratch, radix_bits);

      const double cutoff_report = reached ? cut_cutoff : NA_REAL;

      for (int pi = 0; pi < n_pop; ++pi) {
        const int col = pop_col[pi];
        const int lev = pop_level[pi];

        // Population selection in time order. Overall (col < 0) reuses ord;
        // a subgroup filters ord into sel.
        const int* order_ptr;
        int sz;
        if (col < 0) {
          order_ptr = ord.data();
          sz = m;
        } else {
          int s = 0;
          for (int k = 0; k < m; ++k) {
            const int kk = ord[k];
            if (sub_mat(g0 + orig_cut[kk], col) == lev) sel[s++] = kk;
          }
          order_ptr = sel.data();
          sz = s;
        }

        // Materialize the pooled time-sorted population vectors.
        int n0 = 0, n1 = 0, n_ev = 0, n_drop = 0;
        for (int k = 0; k < sz; ++k) {
          const int kk = order_ptr[k];
          t_sel[k]  = t_cut[kk];
          ei_sel[k] = e_cut[kk];
          j_sel[k]  = j_cut[kk];
          if (j_cut[kk] == 1) ++n1; else ++n0;
          if (e_cut[kk] == 1) ++n_ev;
          if (drop_cut[kk] == 1) ++n_drop;
        }
        const bool both = (n0 > 0 && n1 > 0);

        cutoff_out[pos]  = cutoff_report;
        reached_out[pos] = reached ? 1 : 0;
        nen_out[pos]     = sz;
        nev_out[pos]     = n_ev;
        nd_out[pos]      = n_drop;

        // ---- log-rank (plain / weighted, optionally stratified) ----------
        if (do_logrank && n_ev > 0 && both) {
          double num = NA_REAL, var = NA_REAL;
          if (!use_strata) {
            if (weight_scheme < 0) {
              double r[3];
              logrank_core_impl(t_sel.data(), ei_sel.data(), j_sel.data(), sz, r);
              num = r[0] - r[1];
              var = r[2];
            } else {
              double r[3];
              weighted_logrank_core_impl(t_sel.data(), ei_sel.data(),
                                         j_sel.data(), sz, weight_scheme,
                                         rho, gamma, t_star, r);
              num = r[1];
              var = r[2];
            }
          } else {
            // Reorder the population subset by (stratum, time). order_ptr is
            // already time-sorted; build per-row stratum labels, then sort the
            // index buffer so[] by (stratum, time) with a stable comparator so
            // ties keep time order. Reuse so/ts/es/js/ss as scratch.
            for (int k = 0; k < sz; ++k) {
              st_sel[k] = strata[g0 + orig_cut[order_ptr[k]]];
              so[k] = k;
            }
            std::stable_sort(so.begin(), so.begin() + sz, [&](int p, int q) {
              if (st_sel[p] != st_sel[q]) return st_sel[p] < st_sel[q];
              return t_sel[p] < t_sel[q];
            });
            for (int k = 0; k < sz; ++k) {
              const int kk = so[k];
              ts[k] = t_sel[kk]; es[k] = ei_sel[kk];
              js[k] = j_sel[kk]; ss[k] = st_sel[kk];
            }
            if (weight_scheme < 0) {
              double r[3];
              stratified_logrank_core_impl(ts.data(), es.data(), js.data(),
                                           ss.data(), sz, r);
              num = r[0] - r[1];
              var = r[2];
            } else {
              double r[3];
              stratified_weighted_logrank_core_impl(
                ts.data(), es.data(), js.data(), ss.data(), sz,
                weight_scheme, rho, gamma, t_star, r);
              num = r[1];
              var = r[2];
            }
          }
          lr_mat(pos, 0) = num;
          lr_mat(pos, 1) = var;
        }

        // ---- Cox (PiHE) ---------------------------------------------------
        if (do_coxph && n_ev > 0 && both) {
          double r[4];
          pihe_core_impl(t_sel.data(), ei_sel.data(), j_sel.data(), sz,
                         pihe_ev, r);
          cox_mat(pos, 0) = r[0]; cox_mat(pos, 1) = r[1];
          cox_mat(pos, 2) = r[2]; cox_mat(pos, 3) = r[3];
        }

        // Group splits reused by rmst / km / ahsw / wmst / medsurv.
        const bool need_split = (do_rmst && both) || do_km ||
                                (do_ahsw && n_ev > 0 && both) ||
                                (do_ahr && n_ev > 0 && both) ||
                                (do_wmst && both) ||
                                (do_medsurv && both);
        if (need_split) {
          int a0 = 0, a1 = 0;
          for (int k = 0; k < sz; ++k) {
            if (j_sel[k] == 0) {
              t0d[a0] = t_sel[k]; e0d[a0] = (double) ei_sel[k]; e0i[a0] = ei_sel[k]; ++a0;
            } else {
              t1d[a1] = t_sel[k]; e1d[a1] = (double) ei_sel[k]; e1i[a1] = ei_sel[k]; ++a1;
            }
          }
        }

        // ---- RMST ---------------------------------------------------------
        if (do_rmst && both) {
          double r0[2], r1[2];
          rmst_core_impl(t0d.data(), e0d.data(), n0, tau, rmst_ev, r0);
          rmst_core_impl(t1d.data(), e1d.data(), n1, tau, rmst_ev, r1);
          rmst_mat(pos, 0) = r0[0]; rmst_mat(pos, 1) = r0[1];
          rmst_mat(pos, 2) = r1[0]; rmst_mat(pos, 3) = r1[1];
        }

        // ---- Kaplan-Meier at t_eval (per available group) -----------------
        if (do_km) {
          if (n0 > 0) { double k0[2]; km_core_impl(t0d.data(), e0d.data(), n0, t_eval, k0); km_mat(pos, 0) = k0[0]; }
          if (n1 > 0) { double k1[2]; km_core_impl(t1d.data(), e1d.data(), n1, t_eval, k1); km_mat(pos, 1) = k1[0]; }
        }

        // ---- max-combo (component numerators + covariance) ----------------
        if (do_maxcombo && n_ev > 0 && both) {
          std::fill(mc_ubuf.begin(), mc_ubuf.begin() + nw, 0.0);
          std::fill(mc_vbuf.begin(), mc_vbuf.begin() + (size_t) nw * nw, 0.0);
          combo_logrank_core_impl(t_sel.data(), ei_sel.data(), j_sel.data(), sz,
                                  mc_rho.begin(), mc_gamma.begin(), nw,
                                  mc_ubuf.data(), mc_vbuf.data(), mc_w);
          for (int a = 0; a < nw; ++a) mc_Uout(pos, a) = mc_ubuf[a];
          for (int a = 0; a < nw; ++a)
            for (int b = 0; b < nw; ++b)
              mc_Vout(pos, a * nw + b) = mc_vbuf[(size_t) a * nw + b];
        }

        // ---- AHSW ---------------------------------------------------------
        if (do_ahsw && n_ev > 0 && both) {
          double a0c[6], a1c[6];
          ahsw_core_impl(t0d.data(), e0i.data(), n0, tau,
                         ahsw_dH, ahsw_Rrun, ahsw_Gfrac, a0c);
          ahsw_core_impl(t1d.data(), e1i.data(), n1, tau,
                         ahsw_dH, ahsw_Rrun, ahsw_Gfrac, a1c);
          ahsw_mat(pos, 0) = a0c[2]; ahsw_mat(pos, 1) = a0c[3];
          ahsw_mat(pos, 2) = a0c[4]; ahsw_mat(pos, 3) = (double) n0;
          ahsw_mat(pos, 4) = a1c[2]; ahsw_mat(pos, 5) = a1c[3];
          ahsw_mat(pos, 6) = a1c[4]; ahsw_mat(pos, 7) = (double) n1;
        }

        // ---- Milestone survival (per-group KM + Greenwood at tau) --------
        if (do_milestone && both) {
          double ms[10];
          milestone_core_impl(t_sel.data(), ei_sel.data(), j_sel.data(), sz,
                              tau, ms);
          ms_mat(pos, 0) = ms[0];   // surv0 (control)
          ms_mat(pos, 1) = ms[2];   // var0
          ms_mat(pos, 2) = ms[1];   // surv1 (treatment)
          ms_mat(pos, 3) = ms[3];   // var1
        }

        // ---- Robust modestly-weighted log-rank ---------------------------
        if (do_rmw && n_ev > 0 && both) {
          double rw[6];
          rmw_core_impl(t_sel.data(), ei_sel.data(), j_sel.data(), sz,
                        s_star, rw);
          rmw_mat(pos, 0) = rw[1];  // U_lr
          rmw_mat(pos, 1) = rw[2];  // V_lr
          rmw_mat(pos, 2) = rw[3];  // U_mw
          rmw_mat(pos, 3) = rw[4];  // V_mw
          rmw_mat(pos, 4) = rw[5];  // C
        }

        // ---- Average hazard ratio (Kalbfleisch-Prentice) -----------------
        if (do_ahr && n_ev > 0 && both) {
          double ar[7];
          ahr_core_two_impl(t0d.data(), e0i.data(), n0,
                            t1d.data(), e1i.data(), n1, tau, ar);
          ahr_mat(pos, 0) = ar[0];  // theta1 (control share)
          ahr_mat(pos, 1) = ar[1];  // theta2 (treatment share)
          ahr_mat(pos, 2) = ar[2];  // ahr
          ahr_mat(pos, 3) = ar[3];  // var_theta1
          ahr_mat(pos, 4) = ar[4];  // var_theta2
        }

        // ---- Window mean survival time -----------------------------------
        if (do_wmst && both) {
          double w0[2], w1[2];
          wmst_core_impl(t0d.data(), e0i.data(), n0, wmst_tau1, wmst_tau2,
                         wmst_ev, w0);
          wmst_core_impl(t1d.data(), e1i.data(), n1, wmst_tau1, wmst_tau2,
                         wmst_ev, w1);
          wmst_mat(pos, 0) = w0[0]; wmst_mat(pos, 1) = w0[1];
          wmst_mat(pos, 2) = w1[0]; wmst_mat(pos, 3) = w1[1];
        }

        // ---- Weighted Kaplan-Meier (Pepe-Fleming, pooled) ----------------
        if (do_wkm && n_ev > 0 && both) {
          double wk[4];
          wkm_core_impl(t_sel.data(), ei_sel.data(), j_sel.data(), sz,
                        wkm_weight, wkm_sc, wk);
          wkm_mat(pos, 0) = wk[0];  // num_raw
          wkm_mat(pos, 1) = wk[1];  // variance
          wkm_mat(pos, 2) = wk[2];  // n1 (treatment)
          wkm_mat(pos, 3) = wk[3];  // n2 (control)
        }

        // ---- Median survival (per-group KM + variance summaries) ---------
        if (do_medsurv && both) {
          double bw0 = (medsurv_bw > 0.0) ? medsurv_bw
                                          : silverman_bw(t0d.data(), e0i.data(), n0);
          double bw1 = (medsurv_bw > 0.0) ? medsurv_bw
                                          : silverman_bw(t1d.data(), e1i.data(), n1);
          double m0[6], m1[6];
          medsurv_core_impl(t0d.data(), e0i.data(), n0, bw0, medsurv_sc, m0);
          medsurv_core_impl(t1d.data(), e1i.data(), n1, bw1, medsurv_sc, m1);
          medsurv_mat(pos, 0) = m0[0];  // median (control)
          medsurv_mat(pos, 1) = m0[2];  // greenwood
          medsurv_mat(pos, 2) = m0[3];  // kernel hazard
          medsurv_mat(pos, 3) = m0[4];  // nph variance increment
          medsurv_mat(pos, 4) = m0[5];  // local hazard
          medsurv_mat(pos, 5) = m1[0];  // median (treatment)
          medsurv_mat(pos, 6) = m1[2];
          medsurv_mat(pos, 7) = m1[3];
          medsurv_mat(pos, 8) = m1[4];
          medsurv_mat(pos, 9) = m1[5];
        }

        ++pos;
      }
    }
  }

  return List::create(
    _["cutoff"]     = cutoff_out,
    _["reached"]    = reached_out,
    _["n_enrolled"] = nen_out,
    _["n_event"]    = nev_out,
    _["n_dropout"]  = nd_out,
    _["logrank"]    = do_logrank  ? (SEXP) lr_mat   : R_NilValue,
    _["coxph"]      = do_coxph    ? (SEXP) cox_mat  : R_NilValue,
    _["rmst"]       = do_rmst     ? (SEXP) rmst_mat : R_NilValue,
    _["km"]         = do_km       ? (SEXP) km_mat   : R_NilValue,
    _["ahsw"]       = do_ahsw     ? (SEXP) ahsw_mat : R_NilValue,
    _["mc_U"]       = do_maxcombo ? (SEXP) mc_Uout  : R_NilValue,
    _["mc_V"]       = do_maxcombo ? (SEXP) mc_Vout  : R_NilValue,
    _["milestone"]  = do_milestone ? (SEXP) ms_mat  : R_NilValue,
    _["rmw"]        = do_rmw       ? (SEXP) rmw_mat : R_NilValue,
    _["ahr"]        = do_ahr       ? (SEXP) ahr_mat : R_NilValue,
    _["medsurv"]    = do_medsurv  ? (SEXP) medsurv_mat : R_NilValue,
    _["wkm"]        = do_wkm      ? (SEXP) wkm_mat  : R_NilValue,
    _["wmst"]       = do_wmst     ? (SEXP) wmst_mat : R_NilValue
  );
}
