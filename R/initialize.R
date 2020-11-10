
initialize <- function(data, type_start, hc_init, sigma, omega, theta, gamma, dims)
  # just a wrapper function for mclust::hclass or random initialization
{
  K <- dims$K
  if ( type_start == "hc" ) {
    z <- mclust::unmap( mclust::hclass(hcPairs = hc_init, G = K) )
    # TODO: remove
    # estimate parameters for initial allocation
    #
    # out_mean <- mean_w_array(data, z)
    # inv <- array( diag(q), dim = c(q, q, K) )
    # out_cov_sigma <- cov_w_array(out_mean$data_cent, z, out_mean$mean, inv, out_mean$Nk, pbyp = TRUE)
    # inv <- array( apply(out_cov_sigma$covmat, 3, solve), c(p,p, K) )
    # out_cov_theta <- cov_w_array(out_mean$data_cent, z, out_mean$mean, inv, out_mean$Nk, pbyp = FALSE)

  } else {
    # TODO: random initialization
    tmp <- matrix(runif(dims[[3]]*K), dims[[3]], K)
    z <- sweep(tmp, MARGIN = 1, rowSums(tmp), FUN = "/")
  }
  tau <- colMeans(z)
  ms <- mstep_inverse(data, z, penalty_omega, penalty_gamma, sigma, omega, theta, gamma,
                      control, dims, initialization = TRUE)
  return( list(z = z, parameters = ms$parameters,
               data_cent = ms$data_cent,
               det_sigma = ms$det_sigma, det_theta = ms$det_theta) )
}
