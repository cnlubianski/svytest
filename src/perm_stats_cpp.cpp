// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// Solve weighted least squares: beta = (X'WX)^(-1) X'Wy
inline void wls_fit(const mat& X, const vec& y, const vec& w,
                    vec& beta, vec& mu, double& sigma2) {
  mat WX = X.each_col() % w;         // row-wise multiply columns by w
  mat XtWX = X.t() * WX;
  vec XtWy = X.t() * (w % y);
  beta = solve(XtWX, XtWy);
  mu = X * beta;
  vec resid = y - mu;
  sigma2 = sum(w % square(resid)) / sum(w); // pseudo-ML variance
}

// Compute statistic per type
inline double stat_value(const std::string& stat,
                         const mat& X, const vec& y,
                         const vec& beta, const vec& mu, double sigma2,
                         const vec& w, // current weights
                         const vec& beta0, const vec& mu0, double sigma20,
                         const vec& w_null,
                         const mat& XtX) {
  if (stat == "pred_mean") {
    vec pred = X * beta;
    return mean(pred);
  } else if (stat == "coef_mahal") {
    // Mahalanobis distance with respect to unweighted precision X'X: delta' (X'X) delta
    vec diff = beta - beta0;
    vec tmp  = XtX * diff;
    return dot(diff, tmp);
  } else {
    return NA_REAL;
  }
}

// [[Rcpp::export]]
NumericVector perm_stats_cpp(const arma::mat& X,
                             const arma::vec& y,
                             const arma::vec& w,
                             const arma::vec& w_null,
                             const arma::umat& perm,
                             const std::string& stat) {

  const unsigned int n = X.n_rows;
  const unsigned int B = perm.n_cols;

  // Precompute unweighted XtX (for coef_mahal metric)
  mat XtX = X.t() * X;

  // Baseline equal-weights fit (for coef_dist/coef_mahal/LR reference)
  vec beta0, mu0;
  double sigma20;
  wls_fit(X, y, w_null, beta0, mu0, sigma20);

  NumericVector out(B);

  vec beta, mu;
  double sigma2;

  for (unsigned int b = 0; b < B; ++b) {
    // permuted weights
    vec w_b(n);
    for (unsigned int i = 0; i < n; ++i) {
      w_b[i] = w[ perm(i, b) - 1 ]; // perm is 1-based in R
    }
    // fit
    wls_fit(X, y, w_b, beta, mu, sigma2);
    // stat
    double sb = stat_value(stat, X, y, beta, mu, sigma2, w_b, beta0, mu0, sigma20, w_null, XtX);
    out[b] = sb;
  }

  return out;
}
