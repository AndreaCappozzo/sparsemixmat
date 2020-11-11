#' @export
initialize <- function(data, type_start, hc_init, omega, gamma, dims)
  # just a wrapper function for mclust::hclass or random initialization
{
  K <- dims$K
  if ( type_start == "hc" ) {
    z <- mclust::unmap( mclust::hclass(hcPairs = hc_init, G = K) )
  } else {
    tmp <- matrix(stats::runif(dims[[3]]*K), dims[[3]], K)
    z <- sweep(tmp, MARGIN = 1, rowSums(tmp), FUN = "/")
  }
  tau <- colMeans(z)
  ms <-
    mstep_inverse_sparse_M(
      data = data,
      z = z,
      penalty_omega = penalty_omega,
      penalty_gamma = penalty_gamma,
      penalty_mu = penalty_mu,
      omega = omega,
      gamma = gamma,
      control = control,
      dims = dims,
      initialization = TRUE
    )
  return( list(z = z, parameters = ms$parameters,
               data_cent = ms$data_cent,
               det_sigma = ms$det_sigma, det_theta = ms$det_theta) )
}
