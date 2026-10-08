// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Per-simulation calendar cutoffs for one or more analysis looks.
//
// For look l of simulation s the candidate calendar times are
//   t_cal   = time_looks[l]                            (planned calendar time)
//   t_event = calendar time of the target(s, l)-th counted event
//   t_gap   = cutoff(s, l - 1) + min_gap[l]            (0 for the first look)
//   t_enr   = accrual time of the min_enrolled[l]-th enrolled subject
//             + min_followup[l]
// and the cutoff is
//   min(max(t_cal, t_event, t_gap, t_enr), max_time[l]),
// where an absent condition (NA) is skipped. An event or enrollment target
// that is not met in the simulated data gives +Inf, so the look is not
// reached unless a finite max_time caps it. A look with no lower condition
// is set to max_time[l]. The R wrapper guarantees that every look has at
// least one condition.
//
// The data must be grouped by simulation, with sim_ptr giving 0-based row
// offsets so that simulation s occupies rows [sim_ptr[s], sim_ptr[s + 1]).
// Only rows with count_event == 1 contribute to the event count.
//
// [[Rcpp::export]]
NumericMatrix cutoff_core(const IntegerVector& sim_ptr,
                          const NumericVector& accrual,
                          const NumericVector& tte,
                          const IntegerVector& event,
                          const IntegerVector& count_event,
                          const NumericVector& time_looks,
                          const NumericMatrix& target,
                          const NumericVector& max_time,
                          const NumericVector& min_gap,
                          const NumericVector& min_enrolled,
                          const NumericVector& min_followup) {
  const int nsim    = sim_ptr.size() - 1;
  const int n_looks = time_looks.size();
  NumericMatrix out(nsim, n_looks);

  bool need_enroll = false;
  for (int l = 0; l < n_looks; ++l) {
    if (!std::isnan(min_enrolled[l])) need_enroll = true;
  }

  int max_block = 0;
  for (int s = 0; s < nsim; ++s) {
    const int sz = sim_ptr[s + 1] - sim_ptr[s];
    if (sz > max_block) max_block = sz;
  }
  std::vector<double> cal_ev;  cal_ev.reserve(max_block);
  std::vector<double> acc_sorted; acc_sorted.reserve(max_block);

  for (int s = 0; s < nsim; ++s) {
    const int g0 = sim_ptr[s];
    const int g1 = sim_ptr[s + 1];

    // Calendar times of the counted events, sorted once per simulation.
    cal_ev.clear();
    for (int g = g0; g < g1; ++g) {
      if (event[g] == 1 && count_event[g] == 1) {
        cal_ev.push_back(accrual[g] + tte[g]);
      }
    }
    std::sort(cal_ev.begin(), cal_ev.end());
    const int n_ev = (int) cal_ev.size();

    if (need_enroll) {
      acc_sorted.assign(accrual.begin() + g0, accrual.begin() + g1);
      std::sort(acc_sorted.begin(), acc_sorted.end());
    }
    const int n_acc = g1 - g0;

    double prev = 0.0;
    for (int l = 0; l < n_looks; ++l) {
      bool   any_lower = false;
      double lower     = R_NegInf;

      const double tc = time_looks[l];
      if (!std::isnan(tc)) {
        any_lower = true;
        if (tc > lower) lower = tc;
      }

      const double tg = target(s, l);
      if (!std::isnan(tg)) {
        any_lower = true;
        const int d = (int) tg;
        const double te = (d >= 1 && d <= n_ev) ? cal_ev[d - 1] : R_PosInf;
        if (te > lower) lower = te;
      }

      const double gap = min_gap[l];
      if (!std::isnan(gap)) {
        any_lower = true;
        const double tp = prev + gap;
        if (tp > lower) lower = tp;
      }

      const double ne = min_enrolled[l];
      if (!std::isnan(ne)) {
        any_lower = true;
        const int k = (int) ne;
        const double fu = std::isnan(min_followup[l]) ? 0.0 : min_followup[l];
        const double tn = (k >= 1 && k <= n_acc) ? acc_sorted[k - 1] + fu
                                                 : R_PosInf;
        if (tn > lower) lower = tn;
      }

      double cut = any_lower ? lower : R_PosInf;
      const double cap = max_time[l];
      if (!std::isnan(cap) && cap < cut) cut = cap;

      out(s, l) = cut;
      prev = cut;
    }
  }
  return out;
}
