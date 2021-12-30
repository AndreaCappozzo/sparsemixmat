#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


//' @export
// [[Rcpp::export]]
Rcpp::List penalization_M_mat_lasso(arma::cube data,
                             arma::cube data_cent,
                             arma::colvec z,
                             double Nk,
                             arma::mat mu,
                             arma::mat omega,
                             arma::mat gamma,
                             arma::mat penalty_mu,
                             int p,
                             int q,
                             int N,
                             double CD_tol,
                             int CD_max_iter) {
  
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
  
  bool crit_Q_M =true;
  int iter_Q_M = 0;
  double err_Q_M=0;
  double Q_M =-10000000;
  double Q_M_prev =-10000000;
  
  while(crit_Q_M){
    iter_Q_M += 1;
    
    for(int l=0; l<p; l++){
  // for(int l=(p-1); l>=0; l--){
        arma::uvec ind_no_l = find( row_elem != l );
        arma::uvec ind_l = find( row_elem == l );

      for(int s=0; s<q; s++){
     // for(int s=(q-1); s>=0; s--){
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
    double pen_Q_M = accu(abs(mu)%penalty_mu);
    Q_M = arma::trace(omega*sum_X * gamma * mu.t())- Nk*.5*arma::trace(omega*mu*gamma*mu.t())-pen_Q_M;
    err_Q_M = abs(Q_M - Q_M_prev) / (1 + abs(Q_M));
    Q_M_prev = Q_M;
    crit_Q_M = ((err_Q_M > CD_tol) & (iter_Q_M < CD_max_iter));
  }
  // return mu_penalized;
  return Rcpp::List::create(Named("mu_penalized") = mu_penalized,
                             Named("data_cent_penalized") = data_cent);
  // return Rcpp::List::create(Named("mu_penalized") = mu_penalized,
  //                           Named("Q_M") =Q_M,
  //                            Named("data_cent_penalized") = data_cent);
}


//' @export
// [[Rcpp::export]]
Rcpp::List penalization_M_mat_group_lasso(arma::cube data,
                                           arma::cube data_cent,
                                           arma::colvec z,
                                           double Nk,
                                           arma::mat mu,
                                           arma::mat omega,
                                           arma::mat gamma,
                                           arma::colvec penalty_mu,
                                           int p,
                                           int q,
                                           int N) {
  
  // Containers
  arma::mat mu_penalized=mu;
  arma::mat sum_X(p,q);
  arma::rowvec b_l(q); 
  arma::mat second_addend(p,q);
  arma::mat pinv_omega_lth_col(p,1);
  
  //  Fill containers
  
  sum_X.zeros();
  
  // Vectors of indexes
  arma::colvec row_elem = arma::regspace(0, 1, p-1);
  arma::uvec p_row = find( row_elem <= p );
  
  // weighted sum of X
  for(int i=0; i<N; i++){
    sum_X=sum_X+data.slice(i)*z(i);
  }
  
  arma::mat first_addend = (omega*sum_X)/Nk;
    
  for(int l=0; l<p; l++){

    arma::uvec ind_l = find( row_elem == l );
    
    double norm_lth_row = norm(mu.rows(ind_l), 2);
      // STEP 1: check if the l-th row shall be set to 0
      
      if(norm_lth_row<=penalty_mu(l)){
        mu_penalized.rows(ind_l)*= 0;
      } else {
      // STEP 2: for those rows that are not set to 0, I compute the penalized estimation
  
  second_addend.zeros();
  for(int r=0; r<p; r++){
    
    if (r==l) {
      continue;
    }
    arma::uvec ind_r = find( row_elem == r );
    second_addend+=omega.cols(ind_r)* mu_penalized.rows(ind_r);
  }
  pinv_omega_lth_col = pinv(omega.cols(ind_l).t());
  b_l = pinv_omega_lth_col.t()*(first_addend-second_addend);
  double norm_b_l = norm(b_l, 2);

      // STEP 3: Update rows in mu as per coordinate ascent algorithm
      // mu_penalized.rows(ind_l) =1/(1+penalty_mu(l)/(norm_lth_row-penalty_mu(l)))*b_l;
      mu_penalized.rows(ind_l) =(1-penalty_mu(l)/norm_b_l)*b_l;
      
    }
  }
  // STEP 4: Update data_cent, as per coordinate ascent algorithm
  for(int i=0; i<N; i++){
    data_cent.slice(i)=data.slice(i)-mu_penalized;
  }
  
  // return mu_penalized;
  return Rcpp::List::create(Named("mu_penalized") = mu_penalized,
                            Named("data_cent_penalized") = data_cent);
  
  // return Rcpp::List::create(Named("mu_penalized") = mu_penalized,
  //                           Named("data_cent_penalized") = data_cent,
  //                             Named("first_addend") = first_addend,
  //                             Named("second_addend") = second_addend);
}