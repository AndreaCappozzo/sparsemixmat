#' @export
EM_controls <- function(tol = c(1e-05, sqrt(.Machine$double.eps)),
                        max_iter = rep(1e03, 2),
                        type_start = c("hc", "random"),
                        n_subset_start = NULL,
                        n_random_start = 50,
                        step_width_PGD = 1e-4) # step width for proximal gradient descent (group-lasso penalty only)
  # EM control parameters
{
  list(tol = tol,
       max_iter = max_iter,
       type_start = match.arg( type_start, choices = eval(formals(EM_controls)$type_start) ),
       n_subset_start = n_subset_start,
       n_random_start = n_random_start,
       step_width_PGD=step_width_PGD)
}

#' Scale and centering of three-way objects
#'
#' @param X p  q  N three-way data
#'
#' @return An array of standardized data
#' @export
scale_matrix_data <- function(X, scale=TRUE){

  mean_X <- apply(X, c(1, 2), mean) # p \item q mean matrix
  
  X_cent <- array(apply(X, 3, function(A)
    (A - mean_X)), dim = dim(X))
  
  if (scale == FALSE) {
    return(X_cent)
  }
  
  sd_X <- apply(X, c(1, 2), stats::sd) # p \item q st dev matrix
  
  array(apply(X_cent, 3, function(A)
    (A) / sd_X), dim = dim(X))

}


penalization_M_mat_group_lasso_no_cpp <- function(data = data,
                                                  data_cent,
                                                  z ,
                                                  Nk,
                                                  mu,
                                                  # I initialize it with the sample mean at the current iteration
                                                  omega,
                                                  gamma,
                                                  penalty_mu,
                                                  p = p,
                                                  q = q,
                                                  N = N,
                                                  CD_tol,
                                                  CD_max_iter,
                                                  step_width_PGD) {
  sum_X <- matrix(0, nrow = p, ncol = q)
  mu_penalized <- mu
  
  for (i in 1:N) {
    sum_X <- sum_X + z[i] * data[, , i]
  }
  
  first_addend_gradient <- (omega %*% sum_X%*%gamma)
  
  crit_Q_M <- TRUE
  iter_Q_M <- 0
  Q_M <- Q_M_prev <-  -.Machine$double.xmax
  # sanity check
  Q_M_trace <- c()
  
  while (crit_Q_M) {
    
    iter_Q_M <- iter_Q_M + 1
    
    for (l in 1:p) {
      
      # Proximal gradient method for group lasso with stepsize parameter step_width_PGD
      crit_prox <- TRUE
      mu_pen_l_old <- mu_penalized[l,]
      
      while (crit_prox) {
        
        second_addend_gradient <-
          Nk * (omega %*% mu_penalized %*% gamma) # need to recompute this at each iteration
        gradient_l <-
          -first_addend_gradient[l,] + second_addend_gradient[l,]
        z_l <-
          mu_penalized[l, ] - step_width_PGD * gradient_l # formula 4.16b of \cite{Hastie2015} adapted to our case, with step-size step_width_PGD equal to nu in their formulation
        norm2_z_l <- norm(z_l, type = "2")
        
        if (norm2_z_l < penalty_mu[l]) {
          mu_penalized[l, ] <- 0
        } else{
          mu_penalized[l,] <-
            # (1 - exp(log(penalty_mu[l]) - log(norm2_z_l))) * z_l
            (1 - (step_width_PGD * penalty_mu[l]) / norm2_z_l) * z_l
        }
        mu_pen_l <- mu_penalized[l,]
        crit_prox <- CD_tol<norm(mu_pen_l-mu_pen_l_old,"2")
        mu_pen_l_old <- mu_pen_l
      }
    }
    
    Q_M <-
      sum(diag(omega %*% sum_X %*% gamma %*% t(mu_penalized))) -
      Nk / 2 * sum(diag(omega %*%mu_penalized %*% gamma %*% t(mu_penalized))) -
      sum(penalty_mu * apply(mu_penalized, 1, norm,"2"))
    
    err_Q_M <- abs(Q_M - Q_M_prev) / (1 + abs(Q_M))
    Q_M_prev <- Q_M
    Q_M_trace <- c(Q_M_trace, Q_M)
    
    crit_Q_M <- (err_Q_M > CD_tol & iter_Q_M < CD_max_iter) 
  }
  
  for (i in 1:N) {
    data_cent[,,i] <- data[,,i]-mu_penalized
  }
  # list(mu_penalized=mu_penalized, data_cent_penalized=data_cent)
  list(
    mu_penalized = mu_penalized,
    data_cent_penalized = data_cent,
    Q_M_trace = Q_M_trace,
    Q_M = Q_M
  )
}

penalization_M_mat_coord_ascent_f <- function(type_penalty_mu) {
  switch(type_penalty_mu,
         "lasso" = penalization_M_mat_lasso,
         "group-lasso" = penalization_M_mat_group_lasso)
         # "group-lasso" = penalization_M_mat_group_lasso_no_cpp)
}

pen_mu_f <- function(type_penalty_mu,mu,penalty_mu){
  switch(type_penalty_mu,
         "lasso" = sum( sweep(abs(mu), c(1,2), penalty_mu, "*") ),
         "group-lasso" = sum(sweep(sapply(1:K, function(k)
           apply(mu[, , k], 1, norm, type = "2")),
           MARGIN = 1,
           STATS = penalty_mu,FUN = "*")))
}
