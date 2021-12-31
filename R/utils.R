#' @export
EM_controls <- function(tol = c(1e-05, sqrt(.Machine$double.eps)),
                        max_iter = rep(1e03, 2),
                        type_start = c("hc", "random"),
                        n_subset_start = NULL,
                        n_random_start = 50)
  # EM control parameters
{
  list(tol = tol,
       max_iter = max_iter,
       type_start = match.arg( type_start, choices = eval(formals(EM_controls)$type_start) ),
       n_subset_start = n_subset_start,
       n_random_start = n_random_start)
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
                                                  CD_max_iter) {
  sum_X <- matrix(0, nrow = p, ncol = q)
  mu_penalized <- mu
  
  for (i in 1:N) {
    sum_X <- sum_X + z[i] * data[, , i]
  }
  
  first_addend <- (omega %*% sum_X) / Nk
  ginv_omega <- lapply(1:p, function(l) MASS::ginv(omega[, l])) # Done it once so it needs not be replicated at each iteration in the while loop
  
  crit_Q_M <- TRUE
  iter_Q_M <- 0
  Q_M <- Q_M_prev <-  -.Machine$double.xmax
  # sanity check
  Q_M_trace <- c()
  
  while (crit_Q_M) {
    
    iter_Q_M <- iter_Q_M + 1
    
    for (l in 1:p) {
      second_addend_list <-
        lapply(setdiff(1:p, l), function(r)
          omega[, r] %o% mu_penalized[r,])
      second_addend <- Reduce(f = "+", x = second_addend_list)
      b_l <-
        ginv_omega[[l]] %*% (first_addend - second_addend)
      norm2_b_l <- norm(b_l, type = "2")
      
      if (norm2_b_l < penalty_mu[l]) {
        mu_penalized[l, ] <- 0
      } else{
        mu_penalized[l,] <-
          (1 - exp(log(penalty_mu[l]) - log(norm2_b_l))) * b_l
        # mu_penalized[l, ] <- (1 - penalty_mu[l] /norm(b_l,type = "2")) * b_l
        # (1 - penalty_mu[l] /sqrt(sum(b_l)^2)) * b_l # WRONG
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
