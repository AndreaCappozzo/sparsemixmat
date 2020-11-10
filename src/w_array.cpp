/*
 * Author: Michael Fop
 * Compute weighted means and scattering matrices for array data (N x p x q)
 */
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List mean_w_array(arma::cube data, arma::mat z)
{
  int N = z.n_rows;
  int K = z.n_cols;
  int p = data.n_rows;
  int q = data.n_cols;

  arma::cube mean(p,q,K, fill::zeros);
  arma::vec pro(K, fill::zeros);
  Rcpp::List data_cent(K);

  for ( int k = 0; k < K; k++ ) {
    pro(k) = accu(z.col(k));
    for ( int i = 0; i < N; i++ ) {
      mean.slice(k) += data.slice(i)*z(i,k)/pro(k);
    }
    arma::cube Xc(p,q,N);
    for ( int i = 0; i < N; i++ ) {
      Xc.slice(i) = data.slice(i) - mean.slice(k);
    }
    data_cent[k] = Xc;
  }

  return Rcpp::List::create( Named("Nk") = pro,
                             Named("tau") = pro/N,
                             Named("mean") = mean,
                             Named("data_cent") = data_cent);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List cov_w_array(Rcpp::List data, arma::mat z, arma::cube mean,
                       arma::cube inv, arma::vec Nk,
                       bool pbyp)
{
  int N = z.n_rows;
  int K = z.n_cols;
  int p = mean.n_rows;
  int q = mean.n_cols;

  int c;
  arma::cube covmat;
  if ( pbyp ) {
    covmat.zeros(p,p, K);
    c = q;
  } else {
    covmat.zeros(q,q, K);
    c = p;
  }

  for ( int k = 0; k < K; k++ ) {
    arma::cube data_cent = data[k];
    for ( int i = 0; i < N; i++ ) {
      if ( pbyp ) {
        covmat.slice(k) += data_cent.slice(i) * inv.slice(k) * data_cent.slice(i).t() * z(i,k)/(Nk(k)*c);   // covmat is p x p
      } else {
        covmat.slice(k) += data_cent.slice(i).t() * inv.slice(k) * data_cent.slice(i) * z(i,k)/(Nk(k)*c);  // covmat is q x q
      }
    }
  }

  return Rcpp::List::create( Named("covmat") = covmat );
}