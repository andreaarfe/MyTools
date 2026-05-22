#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <map>
using namespace Rcpp;

// Build surv[k] = P(Binomial(n2, p) >= k), for k = 0..n2+1.
// Entries:
//   surv[0]    = 1.0
//   surv[k]    = R::pbinom(k - 1, n2, p, FALSE, FALSE)  for k = 1..n2
//   surv[n2+1] = 0.0
// Caller can clamp the "need" index into [0, n2+1] to handle boundary cases
// (need <= 0 -> 1.0, need > n2 -> 0.0) without any branching inside hot loops.
static inline void build_surv(int n2, double p, std::vector<double>& surv) {
  surv.assign(n2 + 2, 0.0);
  surv[0] = 1.0;
  for (int k = 1; k <= n2; ++k) {
    surv[k] = R::pbinom(k - 1, n2, p, /*lower_tail=*/0, /*log=*/0);
  }
}

// Compute operating characteristics for one (n1, n, r1, e1, r, p).
// Returns (p_success, en) given precomputed pmf_x1 and surv_x2.
static inline void oc_single(int n1, int n2, int r1, int e1, int r,
                             const std::vector<double>& pmf_x1,
                             const std::vector<double>& surv_x2,
                             double& out_p_success, double& out_en) {
  // P(interim efficacy): X1 in {e1, ..., n1}
  double p_eff1 = 0.0;
  if (e1 <= n1) {
    for (int k = e1; k <= n1; ++k) p_eff1 += pmf_x1[k];
  }

  // P(futility stop): X1 in {0, ..., r1}
  double p_fut1 = 0.0;
  if (r1 >= 0) {
    for (int k = 0; k <= r1; ++k) p_fut1 += pmf_x1[k];
  }

  // Stage-2 success contribution: X1 in {r1+1, ..., e1-1}
  double p_success_2 = 0.0;
  if (r1 + 1 <= e1 - 1) {
    const int hi = n2 + 1;
    for (int x1 = r1 + 1; x1 <= e1 - 1; ++x1) {
      int need = r - x1;
      int idx = need < 0 ? 0 : (need > hi ? hi : need);
      p_success_2 += pmf_x1[x1] * surv_x2[idx];
    }
  }

  out_p_success = p_eff1 + p_success_2;
  double p_stop = p_fut1 + p_eff1;
  out_en = n1 * p_stop + (n1 + n2) * (1.0 - p_stop);
}

// Like oc_single but also returns p_futility and p_efficacy1.
static inline void oc_full(int n1, int n2, int r1, int e1, int r,
                            const std::vector<double>& pmf_x1,
                            const std::vector<double>& surv_x2,
                            double& out_p_success, double& out_p_futility,
                            double& out_p_efficacy1, double& out_en) {
  double p_eff1 = 0.0, p_fut1 = 0.0, p_success_2 = 0.0;
  if (e1 <= n1) {
    for (int k = e1; k <= n1; ++k) p_eff1 += pmf_x1[k];
  }
  if (r1 >= 0) {
    for (int k = 0; k <= r1; ++k) p_fut1 += pmf_x1[k];
  }
  if (r1 + 1 <= e1 - 1) {
    const int hi = n2 + 1;
    for (int x1 = r1 + 1; x1 <= e1 - 1; ++x1) {
      int need = r - x1;
      int idx  = need < 0 ? 0 : (need > hi ? hi : need);
      p_success_2 += pmf_x1[x1] * surv_x2[idx];
    }
  }
  out_p_success   = p_eff1 + p_success_2;
  out_p_futility  = p_fut1;
  out_p_efficacy1 = p_eff1;
  double p_cont   = 1.0 - p_eff1 - p_fut1;
  out_en          = static_cast<double>(n1) + static_cast<double>(n2) * p_cont;
}

// [[Rcpp::export]]
DataFrame find_feasible_designs_cpp(double p0, double p1,
                                    double alpha, double power,
                                    int n_max, int n1_min,
                                    bool irrevocable, bool simon) {
  std::vector<int>    out_n1, out_n, out_n2, out_r1, out_e1, out_r;
  std::vector<double> out_alpha, out_power, out_en0, out_en1;
  std::vector<int>    out_is_irrev;

  // Reusable PMF buffers, sized to the largest n1 we'll see.
  std::vector<double> pmf0(n_max + 1), pmf1(n_max + 1);
  int cached_n1 = -1;

  // Reusable survival buffers for X2 ~ Binomial(n2, p). Rebuilt at each
  // (n, n1) since n2 = n - n1 changes with both. This hoists the R::pbinom
  // calls out of the innermost loop over r/e1/r1/x1.
  std::vector<double> surv0, surv1;

  for (int n = n1_min + 1; n <= n_max; ++n) {
    for (int n1 = n1_min; n1 <= n - 1; ++n1) {
      int n2 = n - n1;

      // Cache PMF of X1 ~ Binomial(n1, p) — independent of r, e1, r1
      if (n1 != cached_n1) {
        for (int k = 0; k <= n1; ++k) {
          pmf0[k] = R::dbinom(k, n1, p0, 0);
          pmf1[k] = R::dbinom(k, n1, p1, 0);
        }
        cached_n1 = n1;
      }

      // Cache survival of X2 ~ Binomial(n2, p) for both p0 and p1.
      build_surv(n2, p0, surv0);
      build_surv(n2, p1, surv1);

      for (int r = 0; r <= n; ++r) {
        if (r == 0) continue;
        if (!simon && irrevocable && r > n1) continue;

        int e1_min = simon ? n1 + 1 : (irrevocable ? r : 0);
        for (int e1 = e1_min; e1 <= n1 + 1; ++e1) {
          for (int r1 = -1; r1 <= e1 - 1; ++r1) {
            double alpha_actual, en_null, power_actual, en_alt;
            oc_single(n1, n2, r1, e1, r, pmf0, surv0, alpha_actual, en_null);
            if (alpha_actual > alpha) continue;
            oc_single(n1, n2, r1, e1, r, pmf1, surv1, power_actual, en_alt);
            if (power_actual < power) continue;

            out_n1.push_back(n1);
            out_n.push_back(n);
            out_n2.push_back(n2);
            out_r1.push_back(r1);
            out_e1.push_back(e1);
            out_r.push_back(r);
            out_alpha.push_back(alpha_actual);
            out_power.push_back(power_actual);
            out_en0.push_back(en_null);
            out_en1.push_back(en_alt);
            out_is_irrev.push_back(e1 >= r ? 1 : 0);
          }
        }
      }
    }
  }

  int m = (int)out_n1.size();
  if (m == 0) {
    return DataFrame::create();
  }

  // Build p0/p1 columns (constant)
  NumericVector p0v(m, p0), p1v(m, p1);
  LogicalVector is_irrev(m);
  for (int i = 0; i < m; ++i) is_irrev[i] = (bool)out_is_irrev[i];

  return DataFrame::create(
    _["n1"]             = out_n1,
    _["n"]              = out_n,
    _["n2"]             = out_n2,
    _["r1"]             = out_r1,
    _["e1"]             = out_e1,
    _["r"]              = out_r,
    _["p0"]             = p0v,
    _["p1"]             = p1v,
    _["alpha_actual"]   = out_alpha,
    _["power_actual"]   = out_power,
    _["en_null"]        = out_en0,
    _["en_alt"]         = out_en1,
    _["is_irrevocable"] = is_irrev,
    _["stringsAsFactors"] = false
  );
}

// 3D Pareto check on (n, en_null, en_alt). Returns logical keep-mask.
//
// Sort-then-sweep: sort by (n, en_null, en_alt) ascending and walk the points
// in lex order, querying a 2D staircase of (en_null, en_alt) over previously
// seen points. Rows with identical (n, en_null, en_alt) form a run that all
// survive together (none dominates the others, since strict-in-at-least-one
// requires a difference). Total work O(m log m).
// [[Rcpp::export]]
LogicalVector find_admissible_designs_cpp(IntegerVector n,
                                          NumericVector en_null,
                                          NumericVector en_alt) {
  int m = n.size();
  LogicalVector keep(m, false);
  if (m == 0) return keep;

  std::vector<int> order(m);
  for (int i = 0; i < m; ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](int a, int b) {
              if (n[a] != n[b]) return n[a] < n[b];
              if (en_null[a] != en_null[b]) return en_null[a] < en_null[b];
              return en_alt[a] < en_alt[b];
            });

  // Staircase: en_null -> min en_alt seen at that key, invariant that as
  // en_null increases along entries, en_alt strictly decreases.
  std::map<double, double> stair;

  int i = 0;
  while (i < m) {
    int j = i;
    int idx0 = order[i];
    int n0 = n[idx0];
    double e0 = en_null[idx0], e1 = en_alt[idx0];
    while (j < m
           && n[order[j]] == n0
           && en_null[order[j]] == e0
           && en_alt[order[j]] == e1) {
      ++j;
    }

    // Dominator query: any staircase entry (k, v) with k <= e0 AND v <= e1?
    bool dominated = false;
    auto it = stair.upper_bound(e0);
    if (it != stair.begin()) {
      --it;
      if (it->second <= e1) dominated = true;
    }

    if (!dominated) {
      for (int k = i; k < j; ++k) keep[order[k]] = true;

      // Insert (e0, e1): remove any (k, v) with k >= e0 AND v >= e1
      // (they are 2D-dominated by the new entry and cannot help future
      // queries that this entry wouldn't already satisfy).
      auto lb = stair.lower_bound(e0);
      while (lb != stair.end() && lb->second >= e1) {
        lb = stair.erase(lb);
      }
      stair[e0] = e1;
    }

    i = j;
  }
  return keep;
}

// [[Rcpp::export]]
DataFrame evaluate_design_cpp(int n1, int n, int r1, int e1, int r,
                              NumericVector p_vec) {
  int m  = p_vec.size();
  int n2 = n - n1;

  NumericVector out_p_success(m), out_p_futility(m),
                out_p_efficacy1(m), out_en(m);

  std::vector<double> pmf(n1 + 1), surv(n2 + 2);

  for (int i = 0; i < m; ++i) {
    double p = p_vec[i];
    for (int x = 0; x <= n1; ++x)
      pmf[x] = R::dbinom(x, n1, p, 0);
    build_surv(n2, p, surv);
    oc_full(n1, n2, r1, e1, r, pmf, surv,
            out_p_success[i], out_p_futility[i],
            out_p_efficacy1[i], out_en[i]);
  }

  return DataFrame::create(Named("p")           = p_vec,
                           Named("p_success")   = out_p_success,
                           Named("p_futility")  = out_p_futility,
                           Named("p_efficacy1") = out_p_efficacy1,
                           Named("en")          = out_en);
}
