#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Per-event-time summary for the WMST computation. Defined at file scope so the
// pointer-based implementation can take a reusable buffer of it from the caller
// and avoid reallocating per cell in a fused simulation loop.
struct WmstEvSummary {
  double e_time;
  double at_risk;
  int n_event;
  double surv_after;
  double tail;
};

// Forward declaration of the pointer-based implementation (defined below).
void wmst_core_impl(const double*, const int*, int, double, double,
                    std::vector<WmstEvSummary>&, double*);

// Single-scan core for the window mean survival time (WMST).
//
// For each group g in 0, ..., ngroup - 1 the function walks the time-sorted
// data once, builds the Kaplan-Meier estimate up to tau2, and returns:
//   column 0: WMST over the window [tau1, tau2], the area under the
//             Kaplan-Meier curve between tau1 and tau2
//   column 1: variance of the WMST estimate
//
// The variance uses the counting process (Greenwood) form
//   var(WMST) = sum_{e_k <= tau2} RW(e_k)^2 d(e_k) / (Y(e_k) (Y(e_k) - d(e_k))),
// where RW(e_k) is the remaining window area integral_{max(e_k, tau1)}^{tau2}
// S(u) du. For tau1 = 0 this reduces to the restricted mean survival time and
// its variance (matching survRM2).
//
// Inputs time, event and grp must be sorted by time in ascending order; grp is
// coded as 0, ..., ngroup - 1; event is 1 for an event and 0 for censoring.

// [[Rcpp::export]]
NumericMatrix wmst_core(NumericVector time, IntegerVector event, IntegerVector grp,
                        int ngroup, double tau1, double tau2) {
  int n = time.size();
  NumericMatrix out(ngroup, 2);

  std::vector<WmstEvSummary> ev;
  std::vector<double> tg;
  std::vector<int> eg;
  tg.reserve(n);
  eg.reserve(n);

  for (int g = 0; g < ngroup; ++g) {
    tg.clear();
    eg.clear();
    for (int i = 0; i < n; ++i) {
      if (grp[i] == g) {
        tg.push_back(time[i]);
        eg.push_back(event[i]);
      }
    }
    double res[2];
    wmst_core_impl(tg.data(), eg.data(), (int) tg.size(), tau1, tau2, ev, res);
    out(g, 0) = res[0];
    out(g, 1) = res[1];
  }

  return out;
}

// Pointer-based implementation with external linkage. The input must be a
// single group, sorted by time in ascending order. Writes wmst, var into
// out[0..1]. The per-event-time summaries are stored in the caller-supplied
// buffer ev, which is cleared on entry and reused across cells. Algorithm
// identical to the per-group body of the exported wrapper above.
void wmst_core_impl(const double* time, const int* event, int n,
                    double tau1, double tau2,
                    std::vector<WmstEvSummary>& ev, double* out) {
  ev.clear();
  if ((int) ev.capacity() < n) ev.reserve(n);

  int atrisk = n;
  double surv = 1.0;

  int i = 0;
  while (i < n) {
    double t0 = time[i];
    if (t0 > tau2) break;
    int d = 0, n_record = 0, j = i;
    while (j < n && time[j] == t0) {
      ++n_record;
      if (event[j] == 1) ++d;
      ++j;
    }
    int y = atrisk;
    if (d > 0 && y > 0) {
      surv *= (1.0 - (double) d / y);
      ev.push_back(WmstEvSummary{t0, (double) y, d, surv, 0.0});
    }
    atrisk -= n_record;
    i = j;
  }

  int m = (int) ev.size();

  double wmst = 0.0;
  double prev_t = 0.0, prev_s = 1.0;
  for (int k = 0; k < m; ++k) {
    double lo = std::max(prev_t, tau1);
    double hi = std::min(ev[k].e_time, tau2);
    if (hi > lo) wmst += prev_s * (hi - lo);
    prev_t = ev[k].e_time;
    prev_s = ev[k].surv_after;
  }
  {
    double lo = std::max(prev_t, tau1);
    double hi = tau2;
    if (hi > lo) wmst += prev_s * (hi - lo);
  }

  if (m >= 1) {
    ev[m - 1].tail = ev[m - 1].surv_after * (tau2 - ev[m - 1].e_time);
    for (int k = m - 2; k >= 0; --k) {
      ev[k].tail = ev[k].surv_after * (ev[k + 1].e_time - ev[k].e_time) +
                   ev[k + 1].tail;
    }
  }

  double variance = 0.0;
  for (int k = 0; k < m; ++k) {
    double rw = (ev[k].e_time >= tau1) ? ev[k].tail : wmst;
    double y = ev[k].at_risk;
    int d = ev[k].n_event;
    if (y - d > 0) {
      variance += rw * rw * (double) d / (y * (y - (double) d));
    }
  }

  out[0] = wmst;
  out[1] = variance;
}
