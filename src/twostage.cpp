#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Compute operating characteristics for one (n1, n, r1, e1, r, p0, p1).
// Helper: returns (p_success, en) for a single response rate p.
static inline void oc_single(int n1, int n2, int r1, int e1, int r,
                             double p, const std::vector<double>& pmf_x1,
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
    for (int x1 = r1 + 1; x1 <= e1 - 1; ++x1) {
      int need = r - x1;
      double p_x2;
      if (need <= 0) {
        p_x2 = 1.0;
      } else if (need > n2) {
        p_x2 = 0.0;
      } else {
        // P(X2 >= need) = pbinom(need - 1, n2, p, lower.tail = FALSE)
        p_x2 = R::pbinom(need - 1, n2, p, /*lower_tail=*/0, /*log=*/0);
      }
      p_success_2 += pmf_x1[x1] * p_x2;
    }
  }

  out_p_success = p_eff1 + p_success_2;
  double p_stop = p_fut1 + p_eff1;
  out_en = n1 * p_stop + (n1 + n2) * (1.0 - p_stop);
}

// [[Rcpp::export]]
List evaluate_design_cpp(int n1, int n, int r1, int e1, int r,
                         double p0, double p1) {
  int n2 = n - n1;

  std::vector<double> pmf0(n1 + 1), pmf1(n1 + 1);
  for (int k = 0; k <= n1; ++k) {
    pmf0[k] = R::dbinom(k, n1, p0, /*log=*/0);
    pmf1[k] = R::dbinom(k, n1, p1, /*log=*/0);
  }

  double alpha_actual, en_null, power_actual, en_alt;
  oc_single(n1, n2, r1, e1, r, p0, pmf0, alpha_actual, en_null);
  oc_single(n1, n2, r1, e1, r, p1, pmf1, power_actual, en_alt);

  return List::create(
    _["n1"]             = n1,
    _["n"]              = n,
    _["n2"]             = n2,
    _["r1"]             = r1,
    _["e1"]             = e1,
    _["r"]              = r,
    _["p0"]             = p0,
    _["p1"]             = p1,
    _["alpha_actual"]   = alpha_actual,
    _["power_actual"]   = power_actual,
    _["en_null"]        = en_null,
    _["en_alt"]         = en_alt,
    _["is_irrevocable"] = (e1 >= r)
  );
}

// [[Rcpp::export]]
DataFrame find_feasible_designs_cpp(double p0, double p1,
                                    double alpha, double power,
                                    int n_max, int n1_min,
                                    bool irrevocable) {
  std::vector<int>    out_n1, out_n, out_n2, out_r1, out_e1, out_r;
  std::vector<double> out_alpha, out_power, out_en0, out_en1;
  std::vector<int>    out_is_irrev;

  // Reusable PMF buffers, sized to the largest n1 we'll see.
  std::vector<double> pmf0(n_max + 1), pmf1(n_max + 1);
  int cached_n1 = -1;

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

      for (int r = 0; r <= n; ++r) {
        if (r == 0) continue;
        if (irrevocable && r > n1) continue;

        int e1_min = irrevocable ? r : 0;
        for (int e1 = e1_min; e1 <= n1 + 1; ++e1) {
          for (int r1 = -1; r1 <= e1 - 1; ++r1) {
            double alpha_actual, en_null, power_actual, en_alt;
            oc_single(n1, n2, r1, e1, r, p0, pmf0, alpha_actual, en_null);
            if (alpha_actual > alpha) continue;
            oc_single(n1, n2, r1, e1, r, p1, pmf1, power_actual, en_alt);
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
// [[Rcpp::export]]
LogicalVector find_admissible_designs_cpp(IntegerVector n,
                                          NumericVector en_null,
                                          NumericVector en_alt) {
  int m = n.size();
  LogicalVector keep(m, true);

  for (int i = 0; i < m; ++i) {
    if (!keep[i]) continue;
    int ni = n[i];
    double e0i = en_null[i], e1i = en_alt[i];
    for (int j = 0; j < m; ++j) {
      if (j == i || !keep[j]) continue;
      if (n[j] <= ni && en_null[j] <= e0i && en_alt[j] <= e1i &&
          (n[j] < ni || en_null[j] < e0i || en_alt[j] < e1i)) {
        keep[i] = false;
        break;
      }
    }
  }
  return keep;
}
