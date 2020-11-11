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
#'
scale_matrix_data <- function(X){

  mean_X <- apply(X, c(1,2), mean) # p \item q mean matrix
  sd_X <- apply(X, c(1,2), stats::sd) # p \item q st dev matrix

  array(apply(X, 3, function(A) (A-mean_X)/sd_X), dim = dim(X))

}
