/*
 * Author: Michael Fop
 * Compute loglikelihood and performs E step of EM algorithm
 */
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List estep_calc(Rcpp::List data, arma::mat z,
                      arma::cube mean, arma::vec tau,
                      arma::cube sigma, arma::cube theta,
                      arma::cube omega, arma::cube gamma,
                      arma::vec det_sigma, arma::vec det_theta)
{
  int N = z.n_rows;
  int K = z.n_cols;
  int p = mean.n_rows;
  int q = mean.n_cols;

  double llk = 0.0;
  double loghood = 0.0;
  double m = 0.0;

  arma::mat dens(N,K), denspro(N,K);

  for ( int k = 0; k < K; k++ ) {
    arma::cube data_cent = data[k];

    for ( int i = 0; i < N; i++ ) {
      dens(i,k) = ( -0.5*trace(
        omega.slice(k) * data_cent.slice(i) * gamma.slice(k) * data_cent.slice(i).t() ) +
          ( -q*0.5*det_sigma(k) -p*0.5*det_theta(k) ) +
          ( -0.5 * p * q * log2pi )
      );
      denspro(i,k) = dens(i,k) + std::log(tau(k));

      if ( k == (K-1) ) {
        m = max( denspro.row(i) );
        loghood = m + std::log( sum( exp(denspro.row(i) - m) ) );
        llk += loghood;
        z.row(i) = exp( denspro.row(i) - loghood );
      }
    }
  }

  return Rcpp::List::create( Named("dens") = dens,
                             Named("denspro") = denspro,
                             Named("z") = z,
                             Named("loglik") = llk );
}
