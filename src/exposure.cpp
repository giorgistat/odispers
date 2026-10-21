#include <Rcpp.h>
using namespace Rcpp;

// Exposure E_i = sum_j len_j w_j k(r_ij / phi_i) / (c phi_i^a) and its derivatives,
// for the Wendland kernel k(t) = (1 - t)^4 (1 + 4t) on [0, 1).
//
// Pairs must be sorted by location and, within each location, by increasing
// distance; start[i] .. start[i + 1] - 1 index the pairs of location i
// (0-based). Because the kernel has compact support, the loop for location i
// stops at the first pair with r >= phi_i, so only pairs within the current
// range are visited.
//
// [[Rcpp::export]]
List exposure_cpp(const IntegerVector& start, const IntegerVector& j,
                  const NumericVector& r, const NumericVector& len,
                  const NumericVector& w, const NumericMatrix& Z,
                  const NumericVector& phi, double a, double c, bool deriv) {
  const int m = phi.size();
  const int q = deriv ? Z.ncol() : 0;
  NumericVector E(m), dE_dlogphi(m);
  NumericMatrix dE_dtheta(m, q);
  for (int i = 0; i < m; ++i) {
    const double ph = phi[i];
    const double scale = 1.0 / (c * std::pow(ph, a));
    double e = 0.0, de = 0.0;
    for (int p = start[i]; p < start[i + 1]; ++p) {
      const double t = r[p] / ph;
      if (t >= 1.0) break;
      const int jj = j[p];
      const double u = 1.0 - t;
      const double u3 = u * u * u;
      const double base = len[jj] * w[jj] * scale;
      const double contrib = base * u3 * u * (1.0 + 4.0 * t);
      e += contrib;
      if (deriv) {
        // d/dlog(phi) of phi^-a k(r/phi) = -phi^-a [a k(t) + t k'(t)]
        de -= base * u3 * (a * u * (1.0 + 4.0 * t) - 20.0 * t * t);
        for (int k = 0; k < q; ++k) dE_dtheta(i, k) += contrib * Z(jj, k);
      }
    }
    E[i] = e;
    dE_dlogphi[i] = de;
  }
  return List::create(_["E"] = E, _["dE_dlogphi"] = dE_dlogphi, _["dE_dtheta"] = dE_dtheta);
}
