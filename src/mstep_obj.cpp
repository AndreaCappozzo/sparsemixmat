/*
 * Author: Michael Fop
 * Compute objective functions in M step of EM algorithm
 */
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
Rcpp::List mstep_obj(Rcpp::List data, arma::mat z, arma::cube mean,
                     arma::cube sigma, arma::cube theta,
                     arma::cube omega, arma::cube gamma, arma::vec tau)
{
  int N = z.n_rows;
  int K = z.n_cols;
  int p = mean.n_rows;
  int q = mean.n_cols;

  double det_sigma, det_theta, sign_sigma, sign_theta, obj;
  arma::vec det_sigma_vec(K), det_theta_vec(K);
  obj = 0.0;

  for ( int k = 0; k < K; k++ ) {
    arma::cube data_cent = data[k];
    log_det(det_sigma, sign_sigma, sigma.slice(k));
    log_det(det_theta, sign_theta, theta.slice(k));
    det_sigma_vec(k) = det_sigma;
    det_theta_vec(k) = det_theta;

    for ( int i = 0; i < N; i++ ) {
      obj += z(i,k) * (log(tau(k)) -0.5*trace(
        omega.slice(k) * data_cent.slice(i) * gamma.slice(k) * data_cent.slice(i).t() ) +
          ( -q*0.5*det_sigma -p*0.5*det_theta )
      );
    }
  }

  return Rcpp::List::create( Named("obj") = obj,
                             Named("det_sigma") = det_sigma_vec,
                             Named("det_theta") = det_theta_vec);
}

