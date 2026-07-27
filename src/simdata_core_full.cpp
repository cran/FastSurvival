// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(dqrng)]]
#include <Rcpp.h>
#include <dqrng.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ------------------------------------------------------------------ //
//  Internal helpers: draw samplers writing into a destination buffer
// ------------------------------------------------------------------ //
// These mirror rpiece_unif_cpp / rpiece_exp_cpp / rcat_cpp exactly, but write
// into a caller-supplied slice so that the whole-pipeline kernel can keep the
// dqrng consumption order identical to the R reference (which calls the
// per-piece samplers in a fixed sequence).

// Piecewise deterministic accrual: each accrual interval receives a fixed
// number of subjects (supplied by the caller as integer counts that sum to the
// group size), repeated every simulation, with uniform entry positions inside
// the interval. One dqrunif draw per subject is consumed for the within-interval
// position only, in subject order, so the later survival and dropout draws see
// the identical dqrng stream and are unaffected.
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

// Exponential / piecewise-exponential survival or dropout times. For a single
// hazard the R reference calls dqrexp(n, rate = hazard); here we draw a unit
// exponential and divide by the rate, which consumes the identical dqrng
// stream and reproduces the same values. For a piecewise hazard this matches
// rpiece_exp_cpp: one unit exponential per subject mapped through the inverse
// cumulative hazard. The piecewise pre-computations (fin_time, cum_haz) are
// supplied by the caller.
static void draw_exp(double* out, int n,
                     const std::vector<double>& hazard,
                     const std::vector<double>& fin_time,
                     const std::vector<double>& cum_haz) {
  const int n_int = (int) hazard.size();
  NumericVector target = dqrng::dqrexp(n, 1.0);
  if (n_int == 1) {
    const double rate = hazard[0];
    for (int i = 0; i < n; ++i) out[i] = target[i] / rate;
    return;
  }
  for (int i = 0; i < n; ++i) {
    const double tgt = target[i];
    int lo = 0, hi = n_int - 1;
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      if (cum_haz[mid] <= tgt) lo = mid + 1; else hi = mid;
    }
    int k = lo;
    if (cum_haz[k] > tgt && k > 0) --k;
    out[i] = fin_time[k] + (tgt - cum_haz[k]) / hazard[k];
  }
}

// Categorical subgroup-cell assignment: one dqrunif draw per subject mapped to
// a 1-based cell via binary search. Matches rcat_cpp.
static void draw_cells(int* out, int n, const std::vector<double>& cum_prev) {
  const int K = (int) cum_prev.size();
  NumericVector u = dqrng::dqrunif(n);
  for (int i = 0; i < n; ++i) {
    const double ui = u[i];
    int lo = 0, hi = K - 1;
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      if (cum_prev[mid] < ui) lo = mid + 1; else hi = mid;
    }
    out[i] = lo + 1;
  }
}

// Convert an R list of numeric vectors to a vector<vector<double>>. Used for
// per-cell hazard / time specifications.
static std::vector<std::vector<double>> to_vv(const List& x) {
  const int m = x.size();
  std::vector<std::vector<double>> out(m);
  for (int i = 0; i < m; ++i) {
    if (x[i] == R_NilValue) { out[i] = std::vector<double>(); continue; }
    NumericVector v = x[i];
    out[i] = std::vector<double>(v.begin(), v.end());
  }
  return out;
}

// ------------------------------------------------------------------ //
//  Simulate one group into pre-allocated column slices
// ------------------------------------------------------------------ //
// Fills the group's block of every output column, starting at row 'base'.
// The dqrng consumption order is: accrual, then (no subgroup) survival,
// dropout; or (subgroup) cell labels, then for each cell in ascending order
// the survival and dropout draws for that cell's subjects in their original
// order. This reproduces the R reference exactly.
//
// Per-cell specs are passed as lists indexed by cell: e_haz[c], e_fin[c],
// e_cum[c] give the survival hazard/pre-computations for cell c, and likewise
// for dropout. For the no-subgroup case n_cell == 1 and cell c = 0 carries the
// single spec.
static void simulate_group_into(
    int base, int nsim, int n, int group_id,
    const std::vector<double>& a_time, const std::vector<int>& acc_counts,
    int n_cell,
    const std::vector<std::vector<double>>& e_haz,
    const std::vector<std::vector<double>>& e_fin,
    const std::vector<std::vector<double>>& e_cum,
    bool has_dropout,
    const std::vector<std::vector<double>>& d_haz,
    const std::vector<std::vector<double>>& d_fin,
    const std::vector<std::vector<double>>& d_cum,
    const std::vector<double>& cum_prev,
    const IntegerMatrix& level_table,   // n_cell x n_fac (0 cols if no subgroup)
    bool fixed_alloc,
    const std::vector<int>& fixed_counts,
    int*    sim_col, int* grp_col,
    double* acc_col, double* srv_col, double* dro_col,
    double* tte_col, int* evt_col, double* cal_col,
    std::vector<int*>& sub_cols         // n_fac pointers (empty if no subgroup)
) {
  const int total_n = nsim * n;

  // sim and group columns.
  for (int s = 0; s < nsim; ++s)
    for (int k = 0; k < n; ++k)
      sim_col[base + s * n + k] = s + 1;
  for (int i = 0; i < total_n; ++i) grp_col[base + i] = group_id;

  // Accrual times for the whole group block (deterministic per-interval counts).
  draw_accrual_det(acc_col + base, nsim, n, a_time, acc_counts);

  if (n_cell == 1) {
    // No subgroups: survival then dropout over the whole block.
    draw_exp(srv_col + base, total_n, e_haz[0], e_fin[0], e_cum[0]);
    if (has_dropout) {
      draw_exp(dro_col + base, total_n, d_haz[0], d_fin[0], d_cum[0]);
    } else {
      for (int i = 0; i < total_n; ++i) dro_col[base + i] = R_PosInf;
    }
    for (int i = 0; i < total_n; ++i) {
      const double sv = srv_col[base + i];
      const double dr = dro_col[base + i];
      const double tt = (sv <= dr) ? sv : dr;
      tte_col[base + i] = tt;
      evt_col[base + i] = (sv <= dr) ? 1 : 0;
      cal_col[base + i] = acc_col[base + i] + tt;
    }
    return;
  }

  // Subgroup path: assign cells, then draw survival/dropout per cell in the
  // subjects' original order, matching the R which(cell == s) loop.
  std::vector<int> cell(total_n);
  if (fixed_alloc) {
    // Deterministic per-cell counts, repeated each simulation; consumes no RNG.
    int pos = 0;
    for (int s = 0; s < nsim; ++s) {
      for (int c = 0; c < n_cell; ++c)
        for (int r = 0; r < fixed_counts[c]; ++r)
          cell[pos++] = c + 1;
    }
  } else {
    draw_cells(cell.data(), total_n, cum_prev);
  }

  // Initialise dropout to Inf (overwritten per cell when dropout is present).
  for (int i = 0; i < total_n; ++i) dro_col[base + i] = R_PosInf;

  // For each cell, gather its subjects (original order), draw survival and
  // dropout for that contiguous count, and scatter back.
  std::vector<int> idx;
  idx.reserve(total_n);
  std::vector<double> buf;
  for (int c = 1; c <= n_cell; ++c) {
    idx.clear();
    for (int i = 0; i < total_n; ++i) if (cell[i] == c) idx.push_back(i);
    const int m = (int) idx.size();
    if (m == 0) continue;

    buf.assign(m, 0.0);
    draw_exp(buf.data(), m, e_haz[c - 1], e_fin[c - 1], e_cum[c - 1]);
    for (int r = 0; r < m; ++r) srv_col[base + idx[r]] = buf[r];

    if (has_dropout) {
      buf.assign(m, 0.0);
      draw_exp(buf.data(), m, d_haz[c - 1], d_fin[c - 1], d_cum[c - 1]);
      for (int r = 0; r < m; ++r) dro_col[base + idx[r]] = buf[r];
    }
  }

  // Derived columns.
  for (int i = 0; i < total_n; ++i) {
    const double sv = srv_col[base + i];
    const double dr = dro_col[base + i];
    const double tt = (sv <= dr) ? sv : dr;
    tte_col[base + i] = tt;
    evt_col[base + i] = (sv <= dr) ? 1 : 0;
    cal_col[base + i] = acc_col[base + i] + tt;
  }

  // Expand the cell label into one subgroup column per factor.
  const int n_fac = level_table.ncol();
  for (int f = 0; f < n_fac; ++f) {
    int* dst = sub_cols[f];
    for (int i = 0; i < total_n; ++i) {
      dst[base + i] = level_table(cell[i] - 1, f);
    }
  }
}

// ------------------------------------------------------------------ //
//  Whole-pipeline simulator
// ------------------------------------------------------------------ //
// Generates the full simulated dataset in one call. Group blocks are written
// in (sim, group) interleaved order directly: for each simulation, the control
// rows then the treatment rows, so no rbind + order is needed afterwards. The
// dqrng draws are issued per group in the same order as the R reference, then
// the rows are placed at their interleaved destinations, which preserves the
// values (placement does not consume RNG).
//
// Note on RNG order vs. row order: the R reference simulates an entire group
// block contiguously (all simulations) and then interleaves. To keep the exact
// dqrng stream, this kernel also draws each group's whole block contiguously
// into a temporary group-contiguous buffer, then copies into the interleaved
// output. This reproduces the reference values bit for bit.
//
// Accrual is deterministic: acc_counts_c and acc_counts_t give the per-interval
// subject counts for the control and treatment groups (each summing to its group
// size); the treatment vector is empty for a one-group simulation.
//
// [[Rcpp::export]]
DataFrame simdata_core_full(
    int nsim,
    const IntegerVector& n_grp,          // length 1 or 2
    const NumericVector& a_time,
    const IntegerVector& acc_counts_c,   // per-interval accrual counts, control
    const IntegerVector& acc_counts_t,   // per-interval accrual counts, treatment
    int n_cell,                          // 1 if no subgroup
    const List& e_haz_c, const List& e_fin_c, const List& e_cum_c,
    const List& e_haz_t, const List& e_fin_t, const List& e_cum_t,
    bool has_dropout,
    const List& d_haz_c, const List& d_fin_c, const List& d_cum_c,
    const List& d_haz_t, const List& d_fin_t, const List& d_cum_t,
    const NumericVector& cum_prev_c, const NumericVector& cum_prev_t,
    const IntegerMatrix& level_table_c, const IntegerMatrix& level_table_t,
    const CharacterVector& sub_names,
    bool fixed_alloc,
    const IntegerVector& fixed_counts_c, const IntegerVector& fixed_counts_t
) {
  const int n_groups = n_grp.size();
  const std::vector<double> at(a_time.begin(), a_time.end());
  const std::vector<int> acc_c(acc_counts_c.begin(), acc_counts_c.end());
  const std::vector<int> acc_t(acc_counts_t.begin(), acc_counts_t.end());

  const auto ehc = to_vv(e_haz_c); const auto efc = to_vv(e_fin_c); const auto ecc = to_vv(e_cum_c);
  const auto eht = to_vv(e_haz_t); const auto eft = to_vv(e_fin_t); const auto ect = to_vv(e_cum_t);
  const auto dhc = to_vv(d_haz_c); const auto dfc = to_vv(d_fin_c); const auto dcc = to_vv(d_cum_c);
  const auto dht = to_vv(d_haz_t); const auto dft = to_vv(d_fin_t); const auto dct = to_vv(d_cum_t);

  const std::vector<double> cpc(cum_prev_c.begin(), cum_prev_c.end());
  const std::vector<double> cpt(cum_prev_t.begin(), cum_prev_t.end());
  const std::vector<int> fcc(fixed_counts_c.begin(), fixed_counts_c.end());
  const std::vector<int> fct(fixed_counts_t.begin(), fixed_counts_t.end());

  const int n_fac = (n_cell > 1) ? level_table_c.ncol() : 0;

  if (n_groups == 1) {
    const int n0 = n_grp[0];
    const int total = nsim * n0;

    IntegerVector sim_col(total), grp_col(total), evt_col(total);
    NumericVector acc_col(total), srv_col(total), dro_col(total),
                  tte_col(total), cal_col(total);
    std::vector<IntegerVector> sub_store(n_fac);
    std::vector<int*> sub_ptr(n_fac);
    for (int f = 0; f < n_fac; ++f) {
      sub_store[f] = IntegerVector(total);
      sub_ptr[f]   = INTEGER(sub_store[f]);
    }

    simulate_group_into(
      0, nsim, n0, 1, at, acc_c, n_cell,
      ehc, efc, ecc, has_dropout, dhc, dfc, dcc,
      cpc, level_table_c, fixed_alloc, fcc,
      INTEGER(sim_col), INTEGER(grp_col),
      REAL(acc_col), REAL(srv_col), REAL(dro_col),
      REAL(tte_col), INTEGER(evt_col), REAL(cal_col), sub_ptr);

    List cols;
    cols["sim"]   = sim_col;
    cols["group"] = grp_col;
    for (int f = 0; f < n_fac; ++f)
      cols[std::string(sub_names[f])] = sub_store[f];
    cols["accrual_time"]  = acc_col;
    cols["surv_time"]     = srv_col;
    cols["dropout_time"]  = dro_col;
    cols["tte"]           = tte_col;
    cols["event"]         = evt_col;
    cols["calendar_time"] = cal_col;
    return DataFrame(cols);
  }

  // Two groups. Draw each group into a group-contiguous temporary (to preserve
  // the dqrng stream), then copy into (sim, group) interleaved positions.
  const int nc = n_grp[0], nt = n_grp[1];
  const int per_sim = nc + nt;
  const int total   = nsim * per_sim;
  const int totc = nsim * nc, tott = nsim * nt;

  // Group-contiguous temporaries.
  IntegerVector simC(totc), grpC(totc), evtC(totc);
  NumericVector accC(totc), srvC(totc), droC(totc), tteC(totc), calC(totc);
  IntegerVector simT(tott), grpT(tott), evtT(tott);
  NumericVector accT(tott), srvT(tott), droT(tott), tteT(tott), calT(tott);
  std::vector<IntegerVector> subC(n_fac), subT(n_fac);
  std::vector<int*> subCptr(n_fac), subTptr(n_fac);
  for (int f = 0; f < n_fac; ++f) {
    subC[f] = IntegerVector(totc); subCptr[f] = INTEGER(subC[f]);
    subT[f] = IntegerVector(tott); subTptr[f] = INTEGER(subT[f]);
  }

  simulate_group_into(
    0, nsim, nc, 1, at, acc_c, n_cell,
    ehc, efc, ecc, has_dropout, dhc, dfc, dcc,
    cpc, level_table_c, fixed_alloc, fcc,
    INTEGER(simC), INTEGER(grpC), REAL(accC), REAL(srvC), REAL(droC),
    REAL(tteC), INTEGER(evtC), REAL(calC), subCptr);

  simulate_group_into(
    0, nsim, nt, 2, at, acc_t, n_cell,
    eht, eft, ect, has_dropout, dht, dft, dct,
    cpt, level_table_t, fixed_alloc, fct,
    INTEGER(simT), INTEGER(grpT), REAL(accT), REAL(srvT), REAL(droT),
    REAL(tteT), INTEGER(evtT), REAL(calT), subTptr);

  // Interleaved output.
  IntegerVector sim_col(total), grp_col(total), evt_col(total);
  NumericVector acc_col(total), srv_col(total), dro_col(total),
                tte_col(total), cal_col(total);
  std::vector<IntegerVector> sub_store(n_fac);
  std::vector<int*> sub_ptr(n_fac);
  for (int f = 0; f < n_fac; ++f) {
    sub_store[f] = IntegerVector(total);
    sub_ptr[f]   = INTEGER(sub_store[f]);
  }

  int* p_sim = INTEGER(sim_col); int* p_grp = INTEGER(grp_col);
  int* p_evt = INTEGER(evt_col);
  double* p_acc = REAL(acc_col); double* p_srv = REAL(srv_col);
  double* p_dro = REAL(dro_col); double* p_tte = REAL(tte_col);
  double* p_cal = REAL(cal_col);

  const int* sC = INTEGER(simC); const int* gC = INTEGER(grpC); const int* eC = INTEGER(evtC);
  const double* aC = REAL(accC); const double* svC = REAL(srvC); const double* dC = REAL(droC);
  const double* tC = REAL(tteC); const double* cC = REAL(calC);
  const int* sT = INTEGER(simT); const int* gT = INTEGER(grpT); const int* eT = INTEGER(evtT);
  const double* aT = REAL(accT); const double* svT = REAL(srvT); const double* dT = REAL(droT);
  const double* tT = REAL(tteT); const double* cT = REAL(calT);

  int out = 0;
  for (int s = 0; s < nsim; ++s) {
    const int cb = s * nc;
    for (int k = 0; k < nc; ++k, ++out) {
      p_sim[out] = sC[cb + k]; p_grp[out] = gC[cb + k]; p_evt[out] = eC[cb + k];
      p_acc[out] = aC[cb + k]; p_srv[out] = svC[cb + k]; p_dro[out] = dC[cb + k];
      p_tte[out] = tC[cb + k]; p_cal[out] = cC[cb + k];
      for (int f = 0; f < n_fac; ++f) sub_ptr[f][out] = subCptr[f][cb + k];
    }
    const int tb = s * nt;
    for (int k = 0; k < nt; ++k, ++out) {
      p_sim[out] = sT[tb + k]; p_grp[out] = gT[tb + k]; p_evt[out] = eT[tb + k];
      p_acc[out] = aT[tb + k]; p_srv[out] = svT[tb + k]; p_dro[out] = dT[tb + k];
      p_tte[out] = tT[tb + k]; p_cal[out] = cT[tb + k];
      for (int f = 0; f < n_fac; ++f) sub_ptr[f][out] = subTptr[f][tb + k];
    }
  }

  List cols;
  cols["sim"]   = sim_col;
  cols["group"] = grp_col;
  for (int f = 0; f < n_fac; ++f)
    cols[std::string(sub_names[f])] = sub_store[f];
  cols["accrual_time"]  = acc_col;
  cols["surv_time"]     = srv_col;
  cols["dropout_time"]  = dro_col;
  cols["tte"]           = tte_col;
  cols["event"]         = evt_col;
  cols["calendar_time"] = cal_col;
  return DataFrame(cols);
}
