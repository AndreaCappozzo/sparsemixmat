#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
Rcpp::List penalization_M_mat_coord_ascent(arma::cube data,
                             arma::cube data_cent,
                             arma::colvec z,
                             double Nk,
                             arma::mat mu,
                             arma::mat omega,
                             arma::mat gamma,
                             arma::mat penalty_mu,
                             int p,
                             int q,
                             int N) {
  
  // Containers
  arma::mat mu_penalized=mu;
  arma::mat mu_penalized_check(1,1);
  arma::mat one_obs_check(1,1);
  arma::mat sum_X(p,q);
  arma::mat numerator_pen_mu(1,1);
  arma::mat denominator_pen_mu(1,1);

  // temp 
  // cube data_cent_copy = data_cent;
  
  //  Fill containers
  mu_penalized.zeros();
  numerator_pen_mu.zeros();
  denominator_pen_mu.zeros();
  sum_X.zeros();

  // Vectors of indexes
  arma::colvec col_elem = arma::regspace(0, 1, q-1);
  arma::colvec row_elem = arma::regspace(0, 1, p-1);
  arma::uvec p_row = find( row_elem <= p );
  arma::uvec q_col = find( col_elem <= q );

  // weighted sum of X
  for(int i=0; i<N; i++){
    sum_X=sum_X+data.slice(i)*z(i);
  }

  for(int l=0; l<p; l++){
      arma::uvec ind_no_l = find( row_elem != l );
      arma::uvec ind_l = find( row_elem == l );

    for(int s=0; s<q; s++){
        arma::uvec ind_no_s = find( col_elem != s ) ;
        arma::uvec ind_s = find( col_elem == s );
        
        // STEP 1: check if the entry [l,s] shall be set to 0
        one_obs_check.zeros();
        mu_penalized_check.zeros();
        
        for(int i=0; i<N; i++){
          // FIXME add eq number from paper
            one_obs_check= omega.submat(ind_l,ind_no_l)*(data_cent.slice(i).submat(ind_no_l,q_col)*gamma.submat(q_col,ind_s))+
                    omega.submat(ind_l,ind_l)*(data_cent.slice(i).submat(ind_l,ind_no_s)*gamma.submat(ind_no_s,ind_s))+
                    omega.submat(ind_l,ind_l)*data.slice(i).submat(ind_l,ind_s)*gamma.submat(ind_s,ind_s);
          
            mu_penalized_check = mu_penalized_check + one_obs_check*z(i);
        }
        
        mu_penalized_check=abs(mu_penalized_check);
        
        if(mu_penalized_check(0,0)<=penalty_mu(l,s)){
          mu_penalized(l,s)=0;
        } else {
        // STEP 2: for those entries that are not set to 0, I compute the penalized MLE
          numerator_pen_mu =
            omega.submat(ind_l,p_row)* sum_X * gamma.submat(q_col,ind_s) - Nk *(
                omega.submat(ind_l,p_row) * mu * gamma.submat(q_col,ind_s) -
                (omega.submat(ind_l,ind_l) * mu.submat(ind_l,ind_s)* gamma.submat(ind_s, ind_s))
            );
          denominator_pen_mu =
            Nk * omega.submat(ind_l,ind_l) * gamma.submat(ind_s,ind_s);
          mu_penalized.submat(ind_l,ind_s)=(numerator_pen_mu-penalty_mu(l,s)*sign(numerator_pen_mu))/
            denominator_pen_mu;
      }
        
      // STEP 3: Update entry [l,s] in mu and data_cent, as per coordinate ascent algorithm
        
      mu.submat(ind_l,ind_s)=mu_penalized.submat(ind_l,ind_s);
        
      for(int i=0; i<N; i++){
        data_cent.slice(i).submat(ind_l,ind_s)=data.slice(i).submat(ind_l,ind_s)-mu_penalized.submat(ind_l,ind_s);
      }
      
    }
  }
  // return mu_penalized;
  return Rcpp::List::create(Named("mu_penalized") = mu_penalized,
                             Named("data_cent_penalized") = data_cent);
}
