#' @export
initialize <-
  function(data,
           type_start,
           hc_init,
           omega,
           gamma,
           dims,
           control,
           penalty_omega,
           penalty_gamma,
           penalty_mu,
           penalization_M_mat_coord_ascent
  )
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
      penalization_M_mat_coord_ascent=penalization_M_mat_coord_ascent
    )
  return( list(z = z, parameters = ms$parameters,
               data_cent = ms$data_cent,
               det_sigma = ms$det_sigma, det_psi = ms$det_psi) )
}
