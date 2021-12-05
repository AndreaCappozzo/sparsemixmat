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
                                                  N = N) {
  sum_X <- matrix(0, nrow = p, ncol = q)
  mu_penalized <- mu
  
  for (i in 1:N) {
    sum_X <- sum_X + z[i] * data[, , i]
  }
  
  second_addend <- (omega %*% sum_X) / Nk
  norm_mu_sample_k <-
    apply(mu, 1, function(b)
      sqrt(sum(b ^ 2)))
  for (l in 1:p) {
    if (norm_mu_sample_k[l] < penalty_mu[l]) {
      mu_penalized[l,] <- 0
    } else{
      first_addend_list <-
        lapply(1:p, function(r)
          omega[, r] %o% mu_penalized[r, ])
      first_addend <- Reduce(f = "+", x = first_addend_list[-l])
      b_l <-
        MASS::ginv(omega[, l]) %*% (-first_addend + second_addend)
      mu_penalized[l, ] <-
        (1 - penalty_mu[l] / (sqrt(sum(b_l) ^ 2))) * b_l
    }
  }
  
  for (i in 1:N) {
    data_cent[,,i] <- data[,,i]-mu_penalized
  }
  list(mu_penalized=mu_penalized, data_cent_penalized=data_cent)
}

penalization_M_mat_coord_ascent_f <- function(type_penalty_mu) {
  switch(type_penalty_mu,
         "lasso" = penalization_M_mat_lasso,
         # "group-lasso" = penalization_M_mat_group_lasso)
         "group-lasso" = penalization_M_mat_group_lasso_no_cpp)
}

pen_mu_f <- function(type_penalty_mu,mu,penalty_mu){
  switch(type_penalty_mu,
         "lasso" = sum( sweep(abs(mu), c(1,2), penalty_mu, "*") ),
         "group-lasso" = sum(sweep(sapply(1:K, function(k)
           apply(mu[, , k], 1, norm, type = "2")),
           MARGIN = 1,
           STATS = penalty_mu,FUN = "*")))
}
