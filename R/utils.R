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


penalization_M_mat_coord_ascent_f <- function(type_penalty_mu) {
  switch(type_penalty_mu,
         "lasso" = penalization_M_mat_lasso,
         "group-lasso" = penalization_M_mat_group_lasso)
}

pen_mu_f <- function(type_penalty_mu,mu,penalty_mu){
  switch(type_penalty_mu,
         "lasso" = sum( sweep(abs(mu), c(1,2), penalty_mu, "*") ),
         "group-lasso" = sum(sweep(sapply(1:K, function(k)
           apply(mu[, , k], 1, norm, type = "2")),
           MARGIN = 1,
           STATS = penalty_mu,FUN = "*")))
}
